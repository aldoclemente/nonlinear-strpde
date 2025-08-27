#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
#include "../include/fe_ls_fisher_kpp.h"
#include "../include/kFoldCV.h"
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;

int main(int argc, char *argv[]){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "../data/simulation_2/";
    std::string mesh_dir = "../data/simulation_1/mesh/";
    
    int sd = std::stoi(argv[1]); // 0 -> "0.00", 1 -> "0.05", 3 -> "0.10" 
    std::string sim = std::string(argv[2]); // 0, ..., 29

    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    matrix_t test_locs = read_mtx<double>("../data/simulation_1/test_locs.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);

    FeSpace Vh(unit_square, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);
    ZeroField<2> f;

    matrix_t time_mesh = read_mtx<double>(data_dir + "time_mesh.mtx");
    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    int n_times = time_mesh.rows();
    int n_sim = 30;
    Triangulation<1, 1> T(time_mesh);

    std::vector<std::string> sigma_lab = {"0.00", "0.05", "0.10"};
    
    std::string sigma_dir = data_dir + "sigma_" + sigma_lab[sd] + "/";
    vector_t rmse = vector_t::Zero(1);

    std::string sim_dir = sigma_dir + sim + "/";

    matrix_t locs = read_mtx<double>(sim_dir + "locs.mtx");
    matrix_t obs = read_mtx<double>(sim_dir + "obs.mtx");

    GeoFrame data(unit_square, T);
    auto &l = data.insert_scalar_layer<POINT, POINT>(
                "layer", std::pair{locs, time_locs});
    l.load_vec("y", obs.reshaped() ); 
            
    vector_t IC = read_mtx<double>(sim_dir + "exact.mtx").col(0);

    double mu = 0.01;
    double alpha = 3.0;
    
    SRPDE nonlinear("y~f", data, fe_ls_fisher_kpp(std::pair(mu,alpha), Vh, IC)); 
    LBFGS<Dynamic> lbfgs;
    int n_dofs = nodes.rows();
    vector_t g_init = 0.5*(vector_t::Random(n_dofs*time_locs.size()) + vector_t::Ones(n_dofs*time_locs.size()));
    //nonlinear.fit(1.0, g_init, lbfgs, BacktrackingLinesearch());
    
    // 10-folds CV to select lambda -----------------------------------------------------------------------------------
    int K = 10;
    KFoldCV cv_(locs, obs, K);
    vector_t cv_error = vector_t::Ones(K);

    ScalarField<1> SSE;

    SSE = [&unit_square, &T, &cv_, &time_locs, &IC, &Vh, &mu, &alpha, &n_times, &n_dofs, &g_init](vector_t lambda) -> double {
        if (lambda[0] <= 0.0) return std::numeric_limits<double>::max();
        double res = 0;
        LBFGS<Dynamic> lbfgs;
        for(int k=0; k < cv_.k(); ++k){

            auto locs = cv_.X_train(k);
            auto obs = cv_.Y_train(k);

            GeoFrame data(unit_square, T);
            auto &l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{locs, time_locs});
            l.load_vec("y", obs.reshaped() ); 
            
            SRPDE model("y~f", data, fe_ls_fisher_kpp(std::pair(mu,alpha), Vh, IC)); 
            
            //vector_t g_init = 0.5*(vector_t::Random(n_dofs*time_locs.size()) + vector_t::Ones(n_dofs*time_locs.size()));
            model.fit(lambda[0], g_init, lbfgs); // BacktrackingLineSearch()

            auto test_locs = cv_.X_test(k);
            auto test_obs = cv_.Y_test(k);
            auto psi_test = internals::point_basis_eval(Vh, test_locs);

            matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
            for(int t = 0; t<n_times-1; ++t){ 
                test_vals.col(t) = psi_test * model.f().block(n_dofs*t, 0, n_dofs, 1);    
            }

            res += std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());
        }
        //std::cout << "lambda " << lambda << " cv_error = " << 1./cv_.k() * res<<  std::endl;
        return 1./cv_.k() * res;

    };

    int n_lambda = 11;
    vector_t lambda_grid = vector_t::Ones(n_lambda);
    for(int i=0; i<n_lambda;++i) lambda_grid(i,0) = std::pow(10, -1. + 0.15 * i);
    
    GridSearch<1> optimizer;
    optimizer.optimize(SSE, lambda_grid); 
    nonlinear.fit(optimizer.optimum()[0], g_init, lbfgs, BacktrackingLineSearch());
    
    /* NelderMead<1> nelder_mead;
    nelder_mead.optimize(SSE,0.5*vector_t::Ones(1));
    
    std::cout << "x opt: " << nelder_mead.optimum()[0]  << std::endl;
    std::cout << "value opt: " << nelder_mead.value() << std::endl;
    std::cout << "#iters: " << nelder_mead.n_iter() << std::endl;

    nonlinear.fit(nelder_mead.optimum()[0], g_init, lbfgs); // , BacktrackingLineSearch() */
    
    Eigen::saveMarket(nonlinear.f(), sim_dir + "estimate_nonlinear.mtx");
    // ----------------------------------------------------------------------------------------------------------------

    // rmse
    auto psi_test = internals::point_basis_eval(Vh, test_locs);

    matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
    for(int t = 0; t<n_times-1; ++t){ 
        test_vals.col(t) = psi_test * nonlinear.f().block(n_dofs*t, 0, n_dofs, 1);    
    }

    matrix_t test_obs = read_mtx<double>(sim_dir + "test_obs.mtx");
    
    rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()); 
    //std::cout << "rmse nonlinear " << rmse << std::endl;
    Eigen::saveMarket(rmse, sim_dir + "rmse_nonlinear.mtx");
    
    std::string command_str = "chown -R 1000:1000 " + data_dir; 
    system(command_str.c_str());
    return 0;
}