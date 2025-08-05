
#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
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
    matrix_t test_locs = read_mtx<double>(data_dir + "test_locs.mtx");
    
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
    int n_lambda = 50;
    matrix_t lambda_grid = matrix_t::Ones(n_lambda,2);
    for(int i=0; i<n_lambda;++i) lambda_grid(i,0) = std::pow(10, -4.0 + 0.1 * i);

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
    
    auto a = integral(unit_square)(mu * dot(grad(u), grad(v)) - alpha * u * v);
    auto F = integral(unit_square)(f * v);
    SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));

    GridSearch<2> optimizer;
    optimizer.optimize(parabolic.gcv(100, 476813), lambda_grid);

    /* std::cout << optimizer.optimum()[0] << " " << optimizer.optimum()[1] << std::endl;
    for(int i=0; i < n_lambda; ++i) std::cout << optimizer.values()[i] << " ";
    std::cout << std::endl;
 */
    parabolic.fit(optimizer.optimum());
    Eigen::saveMarket(locs, sim_dir + "estimate_parabolic.mtx");

    // rmse
    auto psi_test = internals::point_basis_eval(Vh, test_locs);

    int n_dofs = unit_square.nodes().rows();
    matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
    for(int t = 0; t<n_times-1; ++t){ 
        test_vals.col(t) = psi_test * parabolic.f().block(n_dofs*t, 0, n_dofs, 1);    
    }

    matrix_t test_obs = read_mtx<double>(sim_dir + "test_obs.mtx");
    
    rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()); 
   
    Eigen::saveMarket(rmse, sim_dir + "rmse_parabolic.mtx");
    
    std::string command_str = "chown -R 1000:1000 " + data_dir; 
    system(command_str.c_str());

    return 0;
}