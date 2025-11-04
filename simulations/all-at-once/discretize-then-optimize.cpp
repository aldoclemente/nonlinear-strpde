
#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
#include <unsupported/Eigen/SparseExtra>
#include "../include/fe_ls_fisher_kpp.h"
using namespace fdapde;

int main(int argc, char *argv[]){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>::UnitSquare(31, cache_cells);
    
    FeSpace Vh(unit_square, P1<1>);
    matrix_t locs = unit_square.sample(std::stoi(argv[1]));
    auto psi = internals::point_basis_eval(Vh,locs);
    
    matrix_t test_locs = expand_grid(std::vector<vector_t>{
                                         vector_t::LinSpaced(60, 0.05, 0.95), 
                                         vector_t::LinSpaced(60, 0.05, 0.95)});
    auto psi_test = internals::point_basis_eval(Vh, test_locs);

    // Dirichlet BC
    ScalarField<2, decltype([](const PointT& p) { return 0; })> g;
    auto& dof_handler = Vh.dof_handler();
    dof_handler.set_dirichlet_constraint(/* on = */ BoundaryAll, /* data = */ g);

    int n_dofs = dof_handler.n_dofs();
    
    double mu = 1.0;
    double deltaT = 0.001;//5e-1;
    int n_times = 21;
    double Tf = (n_times-1) * deltaT;
    vector_t time_mesh = vector_t::Zero(n_times);
    for(int t = 0; t < n_times; t++) time_mesh[t] = deltaT*t;
    int m = n_times -1 ;    
 
    // define forcing functional
    SpaceTimeField<2, decltype([](const PointT& p, double t) {
        return std::exp(-t)*std::sin(2.*std::numbers::pi*p[0])*std::sin(2.*std::numbers::pi*p[1]);
    })> f_exact;
 
    matrix_t field = matrix_t::Zero(unit_square.nodes().rows(), m+1);
    for(int i = 0; i < n_dofs; ++i){
        for(int j = 0; j < m+1; ++j){
            field(i,j) = f_exact(unit_square.nodes().row(i), time_mesh[j]);
        }
    }

    vector_t IC = field.col(0);

    matrix_t obs = matrix_t::Zero(locs.rows(), m+1);
    for(int j = 0; j < m+1; ++j){
        obs.col(j) = psi*field.col(j);
    }
    matrix_t test_obs = matrix_t::Zero(test_locs.rows(), m+1);
    for(int j = 0; j < m+1; ++j){
        test_obs.col(j) = psi_test*field.col(j); 
    }

    matrix_t misfit = (mu*4.0*std::numbers::pi*std::numbers::pi - 1.0) * field;
    std::cout << "range misfit: " << misfit.minCoeff() << " " << misfit.maxCoeff() << std::endl;
    // --------------------------------------------------------------------------------
    TrialFunction u(Vh);
    TestFunction v(Vh);
   
    Triangulation<1, 1> T(time_mesh);
    matrix_t response = obs.rightCols(m+1);
    vector_t rmse = vector_t::Zero(1);
    
    matrix_t response_parabolic = obs.rightCols(m);
    matrix_t time_locs = time_mesh.tail(m);
    GeoFrame data_para(unit_square, T);
    auto &l_para = data_para.insert_scalar_layer<POINT, POINT>(
                "layer", std::pair{locs, time_locs});
    l_para.load_vec("y", response_parabolic.reshaped() ); 
    
    auto a = integral(unit_square)(mu*dot(grad(u), grad(v)));
    auto F = integral(unit_square)( ZeroField<2>() * v);
    SRPDE parabolic("y ~ f", data_para, fe_ls_parabolic_mono(std::pair{a, F}, IC));

    // ---------------------------------------------------------------------------------------------
    std::cout << "\t running parabolic DtO (all-at-once) ...." << std::endl;
    
    // discretization
    auto bilinear_form = integral(unit_square)(mu*dot(grad(u), grad(v)));
    auto mass = integral(unit_square)(u*v);
    auto A = bilinear_form.assemble();
    auto M = mass.assemble();

    dof_handler.enforce_constraints(A); 
    
    A.makeCompressed();
    M.makeCompressed();
   //auto psi = internals::point_basis_eval(Vh, locs);
    

    // // BUILT time matrices
    double r = 0.0;
    sparse_matrix_t Im(m+1, m+1);   // m x m identity matrix
    Im.setIdentity();
    sparse_matrix_t Mt = kronecker(Im, M); 
    sparse_matrix_t Bt = Mt;
    
    auto At = kronecker(Im, A);
    auto Psi = kronecker(Im, psi);

    sparse_matrix_t M_half = kronecker(Im, M);
    M_half.leftCols(n_dofs) = 0.5*M_half.leftCols(n_dofs);
    M_half.rightCols(n_dofs) = 0.5*M_half.rightCols(n_dofs);
    
    std::cout <<"Mt "<< Mt.rows() << " " << Mt.cols() << std::endl;
    std::cout <<"At " << At.rows() << " " << At.cols() << std::endl;
    std::cout <<"Psi " << Psi.rows() << " " << Psi.cols() << std::endl;
    // assemble matrix associated with derivation in time L_
    // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
    std::vector<Eigen::Triplet<double>> triplet_list;
    triplet_list.reserve((m+1));
    // start assembly loop
    triplet_list.emplace_back(0, 0, 1.0 / deltaT );
    for (int i = 1; i < m+1; ++i) {
        triplet_list.emplace_back(i, i,  1.0 / deltaT );
        triplet_list.emplace_back(i, i - 1, -1.0 / deltaT );
    }
    sparse_matrix_t L(m+1, m+1);   // m x m identity matrix
    L.setFromTriplets(triplet_list.begin(), triplet_list.end());
    L.makeCompressed();

    auto D = kronecker(L, M);
    
    sparse_matrix_t Zero((m+1)*n_dofs, (m+1)*n_dofs);
    Zero.setZero();
    
    // // inizializzarli a zero significa fargli risolvere 
    // // il parabolico monolitico al primo step ?! penso proprio di sì
    
    // //auto Psi_test = kronecker(Im, psi_test);

    // //matrix_t test_obs = read_mtx<double>(mesh_dir + "test_obs.mtx");
    // //matrix_t test_locs = read_mtx<double>(mesh_dir + "test_locs.mtx");
    // //auto psi_test = internals::point_basis_eval(Vh, test_locs);
    int n = locs.rows();

    auto step = [&](double lambda) -> std::pair<vector_t, vector_t> {
        
        vector_t f0 = vector_t::Zero(n_dofs*(m+1));
        vector_t g0 = vector_t::Zero(n_dofs*(m+1));
        vector_t f  = vector_t::Zero(n_dofs*(m+1));
        vector_t g  = vector_t::Zero(n_dofs*(m+1));
        vector_t df = vector_t::Zero(n_dofs*(m+1));
        vector_t dg = vector_t::Zero(n_dofs*(m+1));
        // Sx = b
        sparsesolver_t invS_;
        vector_t x = vector_t::Zero(3 * n_dofs * (m+1));
        vector_t b = vector_t::Zero(3 * n_dofs * (m+1)); // rhs 
        
        bool stop = false;
        int max_iter = 10;
        int iter = 0;
        double tol = 1e-3;
        
        auto PsiTPsi = Psi.transpose()*Psi;
         
        vector_t u_ = vector_t::Zero(n_dofs * (m+1)); //no forcing
        u_.segment(0, n_dofs) += (1.0 / deltaT*M + A)*IC;
        double loss_old = std::numeric_limits<double>::max();
        double loss_new = 0;

        double sse_old = std::numeric_limits<double>::max();
        double sse_new = 0; 
        std::cout << "test vals" << test_obs.rows() << " " << test_obs.cols() << std::endl;

            
            sparse_matrix_t S_00 = 1./(n * (m+1))*Psi.transpose()*Psi;
            sparse_matrix_t S_10 = (D  + At);
            sparse_matrix_t Mt_dir = Mt;

            // right hand side
            b = vector_t::Zero(3 * n_dofs * (m+1));
            b.head(n_dofs * (m+1)) = 1./(n * (m+1)) * Psi.transpose()*response.reshaped();

            b.tail(n_dofs * (m+1)) = u_;
            
            sparse_matrix_t S_01 = S_10.transpose();
            
            auto dirichlet_dofs_ = dof_handler.dirichlet_dofs();

            auto enforce_lhs_dirichlet_bc_ = [&] (SparseBlockMatrix<double, 3, 3>& A) -> void {
                if (dirichlet_dofs_.size() == 0) { return; }
                for( size_t t = 0; t < m +1; ++m){
                for (size_t i = 0; i < dirichlet_dofs_.size(); ++i) {
	            // zero out row and column in correspondance of Dirichlet-type dofs
	            A.row(dirichlet_dofs_[i]+t*n_dofs) *= 0;
	            A.col(dirichlet_dofs_[i]+t*n_dofs) *= 0;
	            A.row(2*t*n_dofs + dirichlet_dofs_[i]) *= 0;
	            A.col(2*t*n_dofs + dirichlet_dofs_[i]) *= 0;
	            // set diagonal elements to 1
	            A.coeffRef(dirichlet_dofs_[i], dirichlet_dofs_[i]) = 1;
	            A.coeffRef(2*t*n_dofs + dirichlet_dofs_[i], 2*t*n_dofs + dirichlet_dofs_[i]) = 1;  
                    }
	            }
                return;
            };

            auto enforce_rhs_dirichlet_bc_ = [&] (vector_t& rhs) -> void {
                if (dirichlet_dofs_.size() == 0) { return; }
                for( size_t t = 0; t < m +1; ++m){
                for (size_t i = 0; i < dirichlet_dofs_.size(); ++i) {
	            // zero out row and column in correspondance of Dirichlet-type dofs
	            rhs.row(dirichlet_dofs_[i]+t*n_dofs) *= 0;
	            rhs.row(2*t*n_dofs + dirichlet_dofs_[i]) *= 0;
	            // set diagonal elements to 1
	            //A.coeffRef(dirichlet_dofs_[i], dirichlet_dofs_[i]) = 1;
	            //A.coeffRef(2*t*n_dofs + dirichlet_dofs_[i], 2*t*n_dofs + dirichlet_dofs_[i]) = 1;  
                    }
	            }
                return;
            };

            SparseBlockMatrix<double, 3, 3> S ( S_00,                    Zero,                 S_01,
                                                Zero,    lambda*deltaT*M_half,      -Mt.transpose(),
                                                S_10,                     -Mt_dir,                 Zero);
            
            enforce_lhs_dirichlet_bc_(S);
            enforce_rhs_dirichlet_bc_(b);

            invS_.compute(S);

            // solve
            x = invS_.solve(b);
            
            // update
            f = x.head(n_dofs * (m+1));
            g = x.block(n_dofs*(m+1), 0, n_dofs*(m+1),1); 

            std::cout << (Psi*f - response.reshaped()).squaredNorm() << std::endl;

            matrix_t test_vals = matrix_t::Zero(test_locs.rows(), (m+1)); 
            for(int t = 0; t < m+1; ++t){ 
                test_vals.col(t) = psi_test * f.block(n_dofs*t, 0, n_dofs, 1);
                std::cout << "\trmse t = " << t << " = " << std::sqrt( (test_vals - test_obs).col(t).array().square().mean()) << std::endl;    
            }
            //test_vals.col(0) = psi_test * IC;    
            
            sse_new = (Psi*f - response.reshaped()).array().square().mean();
            std::cout << "sse f: " << sse_new << std::endl;
            std::cout << "rmse f: " << std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()) << std::endl;
            loss_new = (Psi*f - response.reshaped()).squaredNorm() + lambda*(g.transpose()*(Mt*g))(0,0);
            std::cout << "rmse g: " << std::sqrt( (misfit.reshaped() - g).reshaped().array().square().mean()) << std::endl;
            std::cout << "range f: " << f.minCoeff() << " " << f.maxCoeff() << std::endl; 
            std::cout << "range g: " << g.minCoeff() << " " << g.maxCoeff() << std::endl; 
            std::cout << "J: " <<   loss_new  << std::endl;
            //Eigen::saveMarket(f, "DtO.f.mtx");
            //Eigen::saveMarket(g, "DtO.g.mtx");
        
        return std::make_pair(f, g);
    };

    vector_t lambdas = vector_t::Zero(6,1);
    lambdas << 1e-5, 0.5e-5, 1e-4, 0.5e-4, 0.5e-3, 1e-3;
    //vector_t lambdas = vector_t::Zero(1,1);
    //lambdas << 1e-3;
    for(int i=0; i < lambdas.size(); ++i){
        std::cout << "\t lambda = " << lambdas[i] << std::endl; 
        step(lambdas[i]);
        std::cout << "\t\n parabolic" << std::endl; 
    
        parabolic.fit(std::vector<double>{lambdas[i], 1.0});
        vector_t result = vector_t::Zero((m+1) * n_dofs);
        result.topRows(n_dofs) = IC;
        result.bottomRows(m*n_dofs) = parabolic.f();

        matrix_t test_vals = matrix_t::Zero(test_locs.rows(), (m+1)); 
        for(int t = 0; t < m+1; ++t){ 
            test_vals.col(t) = psi_test * result.block(n_dofs*t, 0, n_dofs, 1);
            std::cout << "\trmse t = " << t << " = " << std::sqrt( (test_vals - test_obs).col(t).array().square().mean()) << std::endl;    
        }

        sparse_matrix_t Im(m,m);
        Im.setIdentity();
        sparse_matrix_t Mt = kronecker(Im,M);
        double sse_new = (Psi*result - response.reshaped()).array().square().mean();
        std::cout << "sse f: " << sse_new << std::endl;
        std::cout << "rmse f: " << std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()) << std::endl;
        std::cout << "rmse g: " << std::sqrt( (misfit.reshaped().head(m*n_dofs) - parabolic.misfit()).reshaped().array().square().mean()) << std::endl;
        std::cout << "range f: " << result.minCoeff() << " " << result.maxCoeff() << std::endl;
        std::cout << "range g: " << parabolic.misfit().minCoeff() << " " << parabolic.misfit().maxCoeff() << std::endl; 
        std::cout << "J: " <<   (Psi*result - response.reshaped()).squaredNorm() + (parabolic.misfit().transpose()*(Mt*parabolic.misfit()))(0,0)  << std::endl; 
        
    }

    Eigen::saveMarket(locs, "locs.mtx");
    Eigen::saveMarket(obs, "obs.mtx");

    return 0;
}