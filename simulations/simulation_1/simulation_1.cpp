
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
    
    std::string data_dir = "input/";
    std::string mesh_dir = "input/mesh/";
    
    std::string n_locs = std::string(argv[1]); // "100", "250" "500" "1000" 
    std::string sim = std::string(argv[2]); // 0, ..., 29
    std::string sim_dir = data_dir + n_locs  + "/" + sim + "/";

    matrix_t locs = read_mtx<double>(sim_dir + "locs.mtx");
    matrix_t obs = read_mtx<double>(sim_dir + "obs_no_noise.mtx");

    //---
    //double min = obs.minCoeff();
    //double range = (obs.maxCoeff() - obs.minCoeff());
    //obs-= min*matrix_t::Ones(obs.rows(), obs.cols());
    //obs/=range;
    obs = obs.unaryExpr([](double x) {return x < 0. ? 0 : x;});
    double min = obs.minCoeff();
    double range = (obs.maxCoeff() - obs.minCoeff());
    std::cout << "min " << min << " range " << range << std::endl;
    //---
    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    matrix_t test_locs = read_mtx<double>(mesh_dir + "test_locs.mtx");
    std::cout << locs.rows() << std::endl;
    std::cout << obs.rows() << " " << obs.cols() << std::endl;
    std::cout << test_locs.rows() << std::endl;
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);
    FeSpace Vh(unit_square, P1<1>);
    //FeSpace Qh(unit_square, P1<1>);
    auto psi = internals::point_basis_eval(Vh,locs);
    

    auto Qh = Vh;
    auto psi_test = internals::point_basis_eval(Vh, test_locs);
    ZeroField<2> f;
    int n_dofs = Vh.dof_handler().n_dofs();
    int n_dofs_control = Qh.dof_handler().n_dofs();
    // IC estimate -------------------------------------------------------------------------------
    TrialFunction u(Vh);
    TestFunction v(Vh);
    
    int n_lambda = 300;
    vector_t exps_IC = vector_t::LinSpaced(n_lambda,-6,0);
    vector_t lambda_grid_IC = vector_t::Ones(n_lambda);
    for(int i=0; i<n_lambda;++i) lambda_grid_IC[i] = std::pow(10, exps_IC[i]);
    
    GeoFrame data_IC(unit_square);
    auto &l_IC = data_IC.insert_scalar_layer<POINT>(
                "layer", locs); 
    l_IC.load_vec("y", obs.col(0).reshaped()); 

    double mu =  1e-3;
    double r = 1.0; 

    auto a_IC = integral(unit_square)(dot(grad(u), grad(v)));
    auto F = integral(unit_square)( ZeroField<2>() * v);
    SRPDE model_IC("y ~ f", data_IC, fe_ls_elliptic(a_IC, F));

    GridSearch<1> optimizer_IC;
    //optimizer_IC.optimize(model_IC.gcv(100, 476813), lambda_grid_IC);
    //model_IC.fit(optimizer_IC.optimum());
    model_IC.fit(1000);
    vector_t IC = model_IC.f(); //exact.col(0);
    matrix_t exact = read_mtx<double>(data_dir + "fisher_kpp_no_lump.mtx");
    std::cout << "IC err: " <<  std::sqrt( ( psi_test*(IC - exact.col(0))).array().square().mean() ) << std::endl; 
    
    // for(int i=0; i<n_dofs; ++i){
    //     if(IC[i] < 0.0) IC[i] = 0;
    // }
    //std::cout << "IC err: " << ( IC - exact.col(0)).array().square().mean() << std::endl; 
    //std::cout << "IC range: " << IC.minCoeff() << " " << IC.maxCoeff()  << std::endl;
    IC = exact.col(0);        
    //std::cout << "IC (exact) " << IC.minCoeff() << " " << IC.maxCoeff()  << std::endl; 
    
    // --------------------------------------------------------------------------------

    matrix_t time_mesh = read_mtx<double>(mesh_dir + "time_mesh.mtx");
    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    double deltaT = time_locs(1,0) - time_locs(0,0);
    int m = time_locs.rows();
    
    Triangulation<1, 1> T(time_mesh);
    matrix_t response = obs.rightCols(m+1);
    vector_t rmse = vector_t::Zero(1);
    GeoFrame data(unit_square, T);
    auto &l = data.insert_scalar_layer<POINT, POINT>(
                "layer", std::pair{locs, time_mesh});
    l.load_vec("y", response.reshaped() ); 
    

    // ---

    matrix_t response_parabolic = obs.rightCols(m);
    GeoFrame data_para(unit_square, T);
    auto &l_para = data_para.insert_scalar_layer<POINT, POINT>(
                "layer", std::pair{locs, time_locs});
    l_para.load_vec("y", response_parabolic.reshaped() ); 
    
    auto a = integral(unit_square)(mu * dot(grad(u), grad(v)) - r * u * v);
    SRPDE parabolic("y ~ f", data_para, fe_ls_parabolic_mono(std::pair{a, F}, IC));

    // matrix_t lambda_grid_time = matrix_t::Ones(n_lambda,2);
    // lambda_grid_time.col(0) = lambda_grid_IC;

    // GridSearch<2> optimizer;
    // //optimizer.optimize(parabolic.gcv(100, 476813), lambda_grid_time);
    // //parabolic.fit(optimizer.optimum());
    
    // vector_t result = vector_t::Zero((m+1) * n_dofs);
    // result.topRows(n_dofs) = IC;
    // //result.bottomRows(m*n_dofs) = parabolic.f();
    
    // std::string out_dir = "output/" + n_locs + "/" + sim + "/";
    // std::string command_str = "mkdir -p " + out_dir; 
    // system(command_str.c_str());

    // Eigen::saveMarket(result, out_dir + "estimate_parabolic.mtx");

    // rmse ---------------------------------------------------------------------

    //auto psi_test = internals::point_basis_eval(Vh, test_locs);
    matrix_t test_vals = matrix_t::Zero(test_locs.rows(), m+1); 

    // for(int t = 0; t < m+1; ++t){ 
    //     test_vals.col(t) = psi_test * result.block(n_dofs*t, 0, n_dofs, 1);    
    // }
    
    matrix_t test_obs = read_mtx<double>(mesh_dir + "test_obs.mtx");
    // std::cout << test_vals.rows() << " " << test_vals.cols() << std::endl;

    // std::cout << test_obs.rows() << " " << test_obs.cols() << std::endl;
    // rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());  
    
    // Eigen::saveMarket(rmse, out_dir + "rmse_parabolic.mtx");

    // std::cout << "rmse para " << rmse << std::endl; 
    
    // --------------------------------------------------------------------------
    //std::cout << "lambda opt (GCV) " << optimizer.optimum()[0] << std::endl;
    SRPDE nonlinear("y~f", data, fe_ls_fisher_kpp(std::pair(mu,r), Vh, Qh, IC,100, 1e-6)); 
    
    //vector_t g_init = read_mtx<double>( mesh_dir + "g_init.mtx").col(2);
    auto& dof_handler_ = Vh.dof_handler();
    
    std::mt19937 gen(12345);
    //vector_t g_init = noise(n_dofs* (m+1), 1, 0.25, gen);
    vector_t g_init = vector_t::Zero(n_dofs_control*(m+1)); // time_locs.size());
    LBFGS<Dynamic> lbfgs(1000, 1e-7, 1e-2, 10);
    //max_iter, tol, step, mem)
    auto SSE2 = [&](vector_t lambda) -> double {
        if (lambda[0] <= 0.0) return std::numeric_limits<double>::max();
        nonlinear.fit(lambda[0], g_init, lbfgs, BacktrackingLineSearch());

        double sse = 1./nonlinear.fitted().rows()*(nonlinear.fitted() - nonlinear.response()).squaredNorm();

        //result.bottomRows(m*n_dofs) = nonlinear.f();
    
        for(int t = 0; t < m+1; ++t){ 
            test_vals.col(t) = psi_test * nonlinear.f().block(n_dofs*t, 0, n_dofs, 1);  
            std::cout << "\trmse t = " << t << " = " << std::sqrt( (test_vals - test_obs).col(t).array().square().mean()) << std::endl;
        }
        rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());

        std::cout << "lambda: " << lambda << " sse: " << sse << " iter# "<< lbfgs.n_iter( ) << " rmse: "<<  rmse << std::endl;
        return sse;
    };
    
    vector_t exps = vector_t::LinSpaced(10,-6,2);
    std::cout << exps << std::endl;
    n_lambda = exps.size();
    vector_t lambda_grid = vector_t::Ones(n_lambda);
    for(int i=0; i<n_lambda;++i) lambda_grid[i] = std::pow(10, exps[i]);
    std::cout << lambda_grid.transpose( ) << std::endl;

    GridSearch<1> optimizer_NL;
    //optimizer_NL.optimize(SSE2, lambda_grid);


    // nonlinear.fit(optimizer.optimum()[0], g_init, lbfgs, BacktrackingLineSearch());
    // std::cout << " n_iter " << lbfgs.n_iter() << std::endl;
    // std::cout<<  nonlinear.misfit().block(n_dofs*(m-1), 0, n_dofs, 1).minCoeff() << " lala " << 
    //              nonlinear.misfit().block(n_dofs*(m-1), 0, n_dofs, 1).maxCoeff() << std::endl;
    // result.bottomRows(m*n_dofs) = nonlinear.f();
    
    // for(int t = 0; t < m+1; ++t){ 
    //     test_vals.col(t) = psi_test * result.block(n_dofs*t, 0, n_dofs, 1);    
    // }
    // rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());  
    // std::cout << "rmse NL " << rmse << std::endl;  
    // Eigen::saveMarket(result, out_dir + "estimate_nonlinear.mtx");
    // Eigen::saveMarket(rmse, out_dir + "rmse_nonlinear.mtx");
    
    // // 
    // SRPDE nonlinear2("y~f", data, fe_ls_fisher_kpp(std::pair(mu,r), Vh, Vh, IC,100, 1e-6)); 
    
    // nonlinear2.fit(optimizer.optimum()[0], parabolic.misfit(), lbfgs, BacktrackingLineSearch());
    
    // result.bottomRows(m*n_dofs) = nonlinear2.f();
    
    // for(int t = 0; t < m+1; ++t){ 
    //     test_vals.col(t) = psi_test * result.block(n_dofs*t, 0, n_dofs, 1);    
    // }
    // rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());  
    // std::cout << "rmse NL2 " << rmse << std::endl; 
    // Eigen::saveMarket(result, out_dir + "estimate_nonlinear_misfit_para.mtx");
    // Eigen::saveMarket(rmse, out_dir + "rmse_nonlinear_misfit_para.mtx");

    // ---------------------------------------------------------------------------------------------
    std::cout << "\t running parabolic DtO (all-at-once) ...." << std::endl;
    
    // discretization
    auto bilinear_form = integral(unit_square)(mu*dot(grad(u), grad(v)));
    auto mass = integral(unit_square)(u*v);
    auto A = bilinear_form.assemble();
    auto M = mass.assemble();
    A.makeCompressed();
    M.makeCompressed();
   //auto psi = internals::point_basis_eval(Vh, locs);
    

    // // BUILT time matrices
    sparse_matrix_t Im(m+1, m+1);   // m x m identity matrix
    Im.setIdentity();
    sparse_matrix_t Mt = kronecker(Im, M); 
    sparse_matrix_t Bt = Mt;
    //Mt.leftCols(n_dofs) = 0.0*M;
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
    std::cout <<"D " << D.rows() << " " << D.cols() << std::endl;
    
    sparse_matrix_t Zero((m+1)*n_dofs, (m+1)*n_dofs);
    Zero.setZero();
    std::cout <<"0 " << Zero.rows() << " " << Zero.cols() << std::endl;
    
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
        u_.segment(0, n_dofs) += (1.0 / deltaT*M + A-r*M)*IC;
        double loss_old = std::numeric_limits<double>::max();
        double loss_new = 0;

        double sse_old = std::numeric_limits<double>::max();
        double sse_new = 0; 
        std::cout << "test vals" << test_obs.rows() << " " << test_obs.cols() << std::endl;

            
            sparse_matrix_t S_00 = 1./(n * (m+1))*Psi.transpose()*Psi;
            sparse_matrix_t S_10 = (D  + At -r*Mt);

            sparse_matrix_t S_01 = S_10.transpose();
            
            SparseBlockMatrix<double, 3, 3> S ( S_00,                    Zero,                 S_01,
                                                Zero,    lambda*deltaT*M_half,      -Bt.transpose(),
                                                S_10,                     -Mt,                 Zero);
            
            
            // Eigen::JacobiSVD<matrix_t> svd(S);
            // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
            // std::cout << "condtion number: " << cond << std::endl;


            invS_.compute(S);

            // right hand side
            b = vector_t::Zero(3 * n_dofs * (m+1));
            b.head(n_dofs * (m+1)) = 1./(n * (m+1)) * Psi.transpose()*response.reshaped();

            b.tail(n_dofs * (m+1)) = u_;
            
            // solve
            x = invS_.solve(b);
            
            // update
            f = x.head(n_dofs * (m+1));
            g = x.block(n_dofs*m, 0, n_dofs*(m+1),1); 

            std::cout << (Psi*f - response.reshaped()).squaredNorm() << std::endl;

            matrix_t test_vals = matrix_t::Zero(test_locs.rows(), (m+1)); 
            for(int t = 0; t < m+1; ++t){ 
                test_vals.col(t) = psi_test * f.block(n_dofs*t, 0, n_dofs, 1);
                std::cout << "\trmse t = " << t << " = " << std::sqrt( (test_vals - test_obs).col(t).array().square().mean()) << std::endl;    
            }
            //test_vals.col(0) = psi_test * IC;    
            
            sse_new = (Psi*f - response.reshaped()).array().square().mean();
            std::cout << "sse: " << sse_new << std::endl;
            std::cout << "rmse: " << std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()) << std::endl;
            loss_new = (Psi*f - response.reshaped()).squaredNorm() + lambda*(g.transpose()*(Mt*g))(0,0);
            std::cout << "range: " << f.minCoeff() << " " << f.maxCoeff() << std::endl; 
        
        
        return std::make_pair(f, g);
    };

    vector_t lambdas = vector_t::Zero(6,1);
    lambdas << 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.;
    for(int i=0; i < lambdas.size(); ++i){
        std::cout << "\t lambda = " << lambdas[i] << std::endl; 
        step(lambdas[i]);
        std::cout << "\t\n parabolic" << std::endl; 
    
        parabolic.fit(std::vector<double>{lambdas[i], 1.0});
        vector_t result = vector_t::Zero((m+1) * n_dofs);
        result.topRows(n_dofs) = IC;
        result.bottomRows(m*n_dofs) = parabolic.f();

        test_vals = matrix_t::Zero(test_locs.rows(), (m+1)); 
        for(int t = 0; t < m+1; ++t){ 
            test_vals.col(t) = psi_test * result.block(n_dofs*t, 0, n_dofs, 1);
            std::cout << "\trmse t = " << t << " = " << std::sqrt( (test_vals - test_obs).col(t).array().square().mean()) << std::endl;    
        }

        double sse_new = (Psi*result - response.reshaped()).array().square().mean();
        std::cout << "sse: " << sse_new << std::endl;
        std::cout << "rmse: " << std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()) << std::endl;
        std::cout << "range: " << result.minCoeff() << " " << result.maxCoeff() << std::endl;
    }

    // ---------------------------------------------------------------------------------------------
//     std::cout << "\t running nonlinear newton based model ...." << std::endl;
    
//     // discretization
//     auto bilinear_form = integral(unit_square)(mu*dot(grad(u), grad(v)));
//     auto mass = integral(unit_square)(u*v);
//     auto A = bilinear_form.assemble();
//     auto M = mass.assemble();
//     A.makeCompressed();
//     M.makeCompressed();
//    //auto psi = internals::point_basis_eval(Vh, locs);
    
//     // // weighted Mass Matrix
//     std::vector<sparse_matrix_t> Mw(m);
//     std::vector<sparse_matrix_t> Gw(m);

//     auto f_old = FeFunction(Vh, vector_t::Zero(n_dofs));
//     auto mass_weight = integral(unit_square)(f_old * u * v);

//     auto update_weighted_mass = [&unit_square, &Vh, &m, &n_dofs, &f_old, &mass_weight](const vector_t& f_k, std::vector<sparse_matrix_t>& Mw) -> void {
//         for(int t=0; t<m+1; ++t){
//             f_old = f_k.block(n_dofs * t, 0, n_dofs, 1);
//             Mw[t] = mass_weight.assemble();
//             Mw[t].makeCompressed();
//         }
//     };

//     // // BUILT time matrices
//     sparse_matrix_t Im(m+1, m+1);   // m x m identity matrix
//     Im.setIdentity();
//     auto Mt = kronecker(Im, M); 
//     auto At = kronecker(Im, A);
//     auto Psi = kronecker(Im, psi);

//     std::cout <<"Mt "<< Mt.rows() << " " << Mt.cols() << std::endl;
//     std::cout <<"At " << At.rows() << " " << At.cols() << std::endl;
//     std::cout <<"Psi " << Psi.rows() << " " << Psi.cols() << std::endl;
//     // assemble matrix associated with derivation in time L_
//     // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
//     std::vector<Eigen::Triplet<double>> triplet_list;
//     triplet_list.reserve(2 * m);
//     // start assembly loop
//     triplet_list.emplace_back(0, 0, 1.0 / deltaT );
//     for (int i = 1; i < m; ++i) {
//         triplet_list.emplace_back(i, i,  1.0 / deltaT );
//         triplet_list.emplace_back(i, i - 1, -1.0 / deltaT );
//     }
//     sparse_matrix_t L(m, m);   // m x m identity matrix
//     L.setFromTriplets(triplet_list.begin(), triplet_list.end());
//     L.makeCompressed();
//     auto D = kronecker(L, M);
    
//     // // inizializzarli a zero significa fargli risolvere 
//     // // il parabolico monolitico al primo step ?! penso proprio di sì
    
//     // //auto Psi_test = kronecker(Im, psi_test);

//     // //matrix_t test_obs = read_mtx<double>(mesh_dir + "test_obs.mtx");
//     // //matrix_t test_locs = read_mtx<double>(mesh_dir + "test_locs.mtx");
//     // //auto psi_test = internals::point_basis_eval(Vh, test_locs);
//     int n = locs.rows();

//     auto step = [&](double lambda) -> std::pair<vector_t, int> {
        
//         vector_t f0 = vector_t::Zero(n_dofs*(m+1));
//         vector_t g0 = vector_t::Zero(n_dofs*(m+1));
//         vector_t f  = vector_t::Zero(n_dofs*(m+1));
//         vector_t g  = vector_t::Zero(n_dofs*(m+1));
//         vector_t df = vector_t::Zero(n_dofs*(m+1));
//         vector_t dg = vector_t::Zero(n_dofs*(m+1));
//         // Sx = b
//         sparsesolver_t invS_;
//         vector_t x = vector_t::Zero(2 * n_dofs * (m+1));
//         vector_t b = vector_t::Zero(2 * n_dofs * (m+1)); // rhs 
        
//         bool stop = false;
//         int max_iter = 10;
//         int iter = 0;
//         double tol = 1e-3;
        
//         auto PsiTPsi = Psi.transpose()*Psi;
        
//         update_weighted_mass(f0, Mw); 
//         update_weighted_mass(g0, Gw); 
        
//         vector_t u_ = vector_t::Zero(n_dofs * (m+1)); //no forcing
//         u_.segment(0, n_dofs) += (1.0 / deltaT) * (M*IC);
//         double loss_old = std::numeric_limits<double>::max();
//         double loss_new = 0;

//         double sse_old = std::numeric_limits<double>::max();
//         double sse_new = 0; 
//         std::cout << "test vals" << test_obs.rows() << " " << test_obs.cols() << std::endl;

//         while( !stop && iter < max_iter){
//             std::cout << "\t iter " << iter << std::endl;
//             update_weighted_mass(f0, Mw); 
//             update_weighted_mass(g0, Gw); 

//             sparse_matrix_t A_tilde = blockDiag(Mw);
//     //         //std::cout << A_tilde.rows() << " " << A_tilde.cols() << std::endl;
//     //         //range(A_tilde);

//             sparse_matrix_t G_tilde = blockDiag(Gw);
//     //         //range(G_tilde);
            
//     //         //std::cout << G_tilde.rows() << " " << G_tilde.cols() << std::endl;
//     //         //std::cout << Psi.rows() << " " << Psi.cols() << std::endl;

//             sparse_matrix_t S_00 = 1./(n * (m+1))*Psi.transpose()*Psi + 2*r*lambda*G_tilde;
//             sparse_matrix_t S_10 = (D  + At -r*Mt + 2*r*A_tilde);

//             sparse_matrix_t S_01 = S_10.transpose();

//             SparseBlockMatrix<double, 2, 2> S ( S_00, lambda*S_01,
//                                                 lambda*S_10, -lambda*Mt);
            
            
//             // Eigen::JacobiSVD<matrix_t> svd(S);
//             // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
//             // std::cout << "condtion number: " << cond << std::endl;


//             invS_.compute(S);

//             // right hand side
//             b = vector_t::Zero(2 * n_dofs * (m+1));
//             b.head(n_dofs * (m+1)) = 1./(n * (m+1)) * Psi.transpose()*response.reshaped() - 1./(n * (m+1))* Psi.transpose()*(Psi*f0) - lambda*S_01*g0 ;
//             b.tail(n_dofs * (m+1)) = -S_10*f0 + Mt*g0;
//             b.tail(n_dofs * (m+1)) += u_; 
//             b.tail(n_dofs * (m+1)) *= lambda;

//             // solve
//             x = invS_.solve(b);
//             df = x.head(n_dofs * (m+1));
//             dg = x.tail(n_dofs * (m+1));
//             // update
//             f = f0 + df;
//             g = g0 + dg;

//             f0 = f;
//             g0 = g;
//             std::cout << (Psi*f - response.reshaped()).squaredNorm() << std::endl;

//             matrix_t test_vals = matrix_t::Zero(test_locs.rows(), (m+1)+1); 
//             for(int t = 0; t < (m+1); ++t){ 
//                 test_vals.col(t+1) = psi_test * f.block(n_dofs*t, 0, n_dofs, 1);    
//             }
//             //test_vals.col(0) = psi_test * IC;    
            
//             sse_new = (Psi*f - response.reshaped()).array().square().mean();
//             std::cout << "sse: " << sse_new << std::endl;
//             std::cout << "rmse: " << std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()) << std::endl;
//             loss_new = (Psi*f - response.reshaped()).squaredNorm() + lambda*(g.transpose()*(Mt*g))(0,0);
//             std::cout << "J      " << loss_new << std::endl;
//             std::cout << "range: " << f.minCoeff() << " " << f.maxCoeff() << std::endl; 
    
//     //         // check

//             //stop = std::sqrt( (df.transpose()*(Mt*df))(0,0) + (dg.transpose()*(Mt*dg))(0,0) ) < tol;
//             stop = std::sqrt( (df.transpose()*(Mt*df))(0,0) + (dg.transpose()*(Mt*dg))(0,0) )  < tol; 
//             std::cout << "stop (Manzoni)    " << std::sqrt( (df.transpose()*(Mt*df))(0,0) + (dg.transpose()*(Mt*dg))(0,0) ) << std::endl;
//             std::cout << "stop (Functional) " << std::abs( loss_new - loss_old ) / loss_old  << std::endl; // / loss_old
//             std::cout << "stop (sse)        " << std::abs( sse_new - sse_old )  /sse_old << std::endl; // / sse_old
//             //stop = std::abs( sse_new - sse_old )  < tol; 
//             sse_old = sse_new;
//             loss_old = loss_new;
//             ++iter;
//         }
        
//         return std::make_pair(f, iter);
//     };


//     std::cout << "\t lambda = 1e-1" << std::endl; 
//     step(1e-1);


//     std::cout << "\t lambda = 1" << std::endl; 
//     step(1);

//     std::cout << "\t lambda = 10" << std::endl; 
//     step(10);

    return 0;
}