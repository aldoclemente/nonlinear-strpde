#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
#include "../include/fe_ls_fisher_kpp.h"
#include "../include/kFoldCV.h"
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;
using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
int main(int argc, char *argv[]){
    

    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "../data/simulation_2/";
    std::string mesh_dir = "../data/simulation_1/mesh/";
    
    int sd = std::stoi(argv[1]); // 0 -> "0.00", 1 -> "0.05", 3 -> "0.10" 
    std::string sim = std::string(argv[2]); // 0, ..., 29
    double lambda = std::stod(argv[3]);
    std::cout << "lambda " << lambda << std::endl;

    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    matrix_t test_locs = read_mtx<double>("../data/simulation_1/test_locs.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);

    FeSpace Vh(unit_square, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);
    
    vector_t time_mesh = read_mtx<double>(data_dir + "time_mesh.mtx");
    vector_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    int m = time_locs.rows();
    double deltaT = time_locs[1] - time_locs[0];
    int n_sim = 30;
    Triangulation<1, 1> T(time_mesh);

    std::vector<std::string> sigma_lab = {"0.00", "0.05", "0.10"};
    
    std::string sigma_dir = data_dir + "sigma_" + sigma_lab[sd] + "/";
    vector_t rmse = vector_t::Zero(1);

    std::string sim_dir = sigma_dir + sim + "/";

    matrix_t locs = read_mtx<double>(sim_dir + "locs.mtx");
    vector_t obs = read_mtx<double>(sim_dir + "obs.mtx").reshaped();
    std::cout << "obd " <<  obs.rows() << std::endl;
    int n = locs.rows();
    
    std::cout << time_locs.transpose() << std::endl; 
    //DEVI METTERE LA IC al termine noto

    vector_t IC = read_mtx<double>(sim_dir + "ic.mtx");
    std::cout << IC.rows() << std::endl;
    //vector_t IC = read_mtx<double>(sim_dir + "exact.mtx").col(0);
    double mu = 0.01; // 0.001
    double r = 1.0; // 1.0


    // discretization
    auto bilinear_form = integral(unit_square)(mu*dot(grad(u), grad(v)));
    auto mass = integral(unit_square)(u*v);
    auto A = bilinear_form.assemble();
    auto M = mass.assemble();
    A.makeCompressed();
    M.makeCompressed();
    auto psi = internals::point_basis_eval(Vh, locs);
    
    auto& dof_handler = Vh.dof_handler();
    auto n_dofs = dof_handler.n_dofs();

    // weighted Mass Matrix
    std::vector<sparse_matrix_t> Mw(m);
    std::vector<sparse_matrix_t> Gw(m);

    auto f_old = FeFunction(Vh, vector_t::Zero(n_dofs));
    auto mass_weight = integral(unit_square)(f_old * u * v);

    auto update_weighted_mass = [&unit_square, &Vh, &m, &n_dofs, &f_old, &mass_weight](const vector_t& f_k, std::vector<sparse_matrix_t>& Mw) -> void {
        for(int t=0; t<m; ++t){
            f_old = f_k.block(n_dofs * t, 0, n_dofs, 1);
            Mw[t] = mass_weight.assemble();
            Mw[t].makeCompressed();
        }
    };

    // BUILT time matrices
    sparse_matrix_t Im(m, m);   // m x m identity matrix
    Im.setIdentity();
    auto Mt = kronecker(Im, M); 
    auto At = kronecker(Im, A);
    auto Psi = kronecker(Im, psi);

    std::cout <<"Mt "<< Mt.rows() << " " << Mt.cols() << std::endl;
    std::cout <<"At " << At.rows() << " " << At.cols() << std::endl;
    std::cout <<"Psi " << Psi.rows() << " " << Psi.cols() << std::endl;
    // assemble matrix associated with derivation in time L_
    // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
    std::vector<Eigen::Triplet<double>> triplet_list;
    triplet_list.reserve(2 * m);
    // start assembly loop
    triplet_list.emplace_back(0, 0, 1.0 / deltaT );
    for (int i = 1; i < m; ++i) {
        triplet_list.emplace_back(i, i,  1.0 / deltaT );
        triplet_list.emplace_back(i, i - 1, -1.0 / deltaT );
    }
    sparse_matrix_t L(m, m);   // m x m identity matrix
    L.setFromTriplets(triplet_list.begin(), triplet_list.end());
    L.makeCompressed();
    auto D = kronecker(L, M);
    
    // inizializzarli a zero significa fargli risolvere 
    // il parabolico monolitico al primo step ?! penso proprio di sì
    
    vector_t f0 = vector_t::Zero(n_dofs*m);
    vector_t g0 = vector_t::Zero(n_dofs*m);
    vector_t f  = vector_t::Zero(n_dofs*m);
    vector_t g  = vector_t::Zero(n_dofs*m);
    vector_t df = vector_t::Zero(n_dofs*m);
    vector_t dg = vector_t::Zero(n_dofs*m);
    // Sx = b
    sparsesolver_t invS_;
    vector_t x = vector_t::Zero(2 * n_dofs * m);
    vector_t b = vector_t::Zero(2 * n_dofs * m); // rhs 
    
    bool stop = false;
    int max_iter = 10;
    int iter = 0;
    double tol = 1e-4;
    
    auto PsiTPsi = Psi.transpose()*Psi;
    std::cout << "PsiTPsi " << PsiTPsi.rows() <<  " " << PsiTPsi.cols() << std::endl;
    update_weighted_mass(f0, Mw); 
    update_weighted_mass(g0, Gw); 
    std::cout << "f0 " << f0.rows() << std::endl;
    std::cout << "Mw " <<  Mw[0].rows() << " " << Mw[0].cols() << std::endl;

    vector_t u_ = vector_t::Zero(n_dofs * m); //no forcing
    u_.segment(0, n_dofs) += (1.0 / deltaT) * (M*IC);

    while( !stop && iter < max_iter){
        std::cout << "\t iter " << iter  << std::endl;
        update_weighted_mass(f0, Mw); 
        update_weighted_mass(g0, Gw); 

        sparse_matrix_t A_tilde = blockDiag(Mw);
        std::cout << A_tilde.rows() << " " << A_tilde.cols() << std::endl;
        range(A_tilde);

        sparse_matrix_t G_tilde = blockDiag(Gw);
        range(G_tilde);
        
        std::cout << G_tilde.rows() << " " << G_tilde.cols() << std::endl;
        std::cout << Psi.rows() << " " << Psi.cols() << std::endl;

        sparse_matrix_t S_00 = 1./(n * m)*Psi.transpose()*Psi + 2*r*lambda*G_tilde;
        sparse_matrix_t S_10 = (D  + At -r*Mt + 2*r*A_tilde);

        sparse_matrix_t S_01 = lambda*S_10.transpose();

        SparseBlockMatrix<double, 2, 2> S ( S_00, S_01,
                                            S_10, -Mt);
        
        
        invS_.compute(S);

        // right hand side
        b = vector_t::Zero(2 * n_dofs * m);
        b.head(n_dofs * m) = 1./(n * m) * Psi.transpose()*obs - 1./(n * m)* Psi.transpose()*(Psi*f0) - lambda*S_01*g0 ;
        b.tail(n_dofs * m) = -S_10*f0 + Mt*g0;
        b.tail(n_dofs * m) += u_; 

         // solve
        x = invS_.solve(b);
        df = x.head(n_dofs * m);
        dg = x.tail(n_dofs * m);
        // update
        f = f0 + df;
        g = g0 + dg;

        std::cout << "error: " << std::sqrt( (Psi*f - obs).array().square().mean() ) << std::endl;
        f0 = f;
        g0 = g;

        // check
        stop = std::sqrt( (df.transpose()*(Mt*df))(0,0) + (dg.transpose()*(Mt*dg))(0,0) ) < tol;
        std::cout << "stop? " << std::sqrt( (df.transpose()*(Mt*df))(0,0) + (dg.transpose()*(Mt*dg))(0,0) ) << std::endl;
        ++iter;
    }

    std::cout << "#iter "<< iter << std::endl;

    auto psi_test = internals::point_basis_eval(Vh, test_locs);
    matrix_t test_obs = read_mtx<double>(sim_dir + "test_obs.mtx");

    matrix_t test_vals = matrix_t::Zero(test_locs.rows(), m); // t0 butto via
    for(int t = 0; t < m; ++t){ 
        test_vals.col(t) = psi_test * f.block(n_dofs*t, 0, n_dofs, 1);    
    }

    rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()); 
    std::cout << "rmse: " << rmse << std::endl;
    
    return 0;

}