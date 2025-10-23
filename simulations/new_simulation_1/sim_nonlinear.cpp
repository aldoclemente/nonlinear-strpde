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

    std::string data_dir = "data";
    std::string mesh_dir = "data/mesh/";
    
    int n_locs = std::stoi(argv[1]); // 0 -> "100", 1 -> "250", 3 -> "500", 4 "1000" 
    std::string sim = std::string(argv[2]); // 0, ..., 29

    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);
    FeSpace Vh(unit_square, P1<1>);   // piecewise linear continuous scalar finite elements
    
    std::string command_str;

    vector_t time_mesh = read_mtx<double>(mesh_dir + "time_mesh_pde_sol.mtx");
    int n_times = time_mesh.size();
    
    std::vector<std::string> sigma_lab = {"100", "250", "500", "1000"};
    std::string sigma_dir = data_dir + "/" + sigma_lab[sd] + "/";
    vector_t rmse = vector_t::Zero(1);
    std::string sim_dir = sigma_dir + sim + "/";

    matrix_t locs = read_mtx<double>(sim_dir + "locs.mtx");
    matrix_t obs = read_mtx<double>(sim_dir + "obs.mtx");

    // IC estimate -------------------------------------------------------------------------------
    int n_lambda = 50
    auto lambda_grid = vector_t::Ones(n_lambda);
    for(int i=0; i<n_lambda;++i) lambda_grid(i,0) = std::pow(10, -4. + 0.15 * i);

    GeoFrame data_IC(unit_square);
    auto &l_IC = data.insert_scalar_layer<POINT>(
                "layer", locs); 
    l.load_vec("y", obs.col(0).reshaped()); 

    double mu =  1.;// 0.01;
    double alpha = 10.0; //3.0;

    auto a = integral(unit_square)(mu*dot(grad(u), grad(v)));
    auto F = integral(unit_square)(f * v);
    SRPDE model_IC("y ~ f", data, fe_ls_elliptic(a, F));

    GridSearch<1> optimizer_IC;
    optimizer_IC.optimize(model_IC.gcv(100, 476813), lambda_grid);
    model_IC.fit(optimizer.optimum());
    vector_t IC = model_IC.f();
    // -------------------------------------------------------------------------------------------
    
    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    int n_times = time_mesh.rows();
    int n_sim = 30;
    Triangulation<1, 1> T(time_mesh);

    GeoFrame data(unit_square, T);
    auto &l = data.insert_scalar_layer<POINT, POINT>(
                "layer", std::pair{locs, time_locs});
    l.load_vec("y", obs.rightCols(n_times-1).reshaped() ); 
            
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
            l.load_vec("y", obs.rightCols(n_times-1).reshaped() ); 
            
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
    std::cout << "lambda: " << lambda_grid.transpose() << std::endl;
    std::cout <<"opt: "<<  optimizer.optimum() << std::endl;
    std::cout << optimizer.values() << std::endl;
    nonlinear.fit(optimizer.optimum()[0], g_init, lbfgs, BacktrackingLineSearch());
    Eigen::saveMarket( optimizer.optimum(), sim_dir + "lambda_nonlinear.mtx");
    /* NelderMead<1> nelder_mead;
    nelder_mead.optimize(SSE,0.5*vector_t::Ones(1));
    
    std::cout << "x opt: " << nelder_mead.optimum()[0]  << std::endl;
    std::cout << "value opt: " << nelder_mead.value() << std::endl;
    std::cout << "#iters: " << nelder_mead.n_iter() << std::endl;

    nonlinear.fit(nelder_mead.optimum()[0], g_init, lbfgs); // , BacktrackingLineSearch() */
    vector_t result = vector_t::Zero(n_times * n_dofs);
    result.topRows(n_dofs) = IC;
    result.bottomRows((n_times-1)*(n_dofs)) = nonlinear.f();
    Eigen::saveMarket(result, sim_dir + "estimate_nonlinear.mtx");
    // ----------------------------------------------------------------------------------------------------------------

    // rmse
    auto psi_test = internals::point_basis_eval(Vh, test_locs);

    matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
    for(int t = 0; t<n_times; ++t){ 
        test_vals.col(t) = psi_test * result.block(n_dofs*t, 0, n_dofs, 1);    
    }

    matrix_t test_obs = read_mtx<double>(sim_dir + "test_obs.mtx");
    
    rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()); 
    //std::cout << "rmse nonlinear " << rmse << std::endl;
    Eigen::saveMarket(rmse, sim_dir + "rmse_nonlinear.mtx");
    std::cout << "rmse: " << rmse << std::endl;
    std::string command_str = "chown -R 1000:1000 " + data_dir; 
    system(command_str.c_str());
    return 0;
}