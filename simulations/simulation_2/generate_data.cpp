

#include <fdaPDE/core/fdaPDE/core.h>
#include <fdaPDE/src/solvers/utility.h>
#include "../include/utils.h"
#include <unsupported/Eigen/SparseExtra>
#include <string>
#include <vector>

auto noise(std::size_t n, double sigma, std::mt19937 gen) {
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n, 1);
    std::normal_distribution<> __noise(0.0, sigma);
    for (std::size_t i = 0; i < n; ++i) { res(i, 0) = __noise(gen); }
    return res;
}

using namespace fdapde;
int main (){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "../data/simulation_2/";
    std::string mesh_dir = "../data/simulation_1/mesh/";
    
    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);
    FeSpace Vh(unit_square, P1<1>);   // piecewise linear continuous scalar finite elements
 
    // To check if the directory exist or not, create it if
    // doesn't exist
    matrix_t time_mesh = read_mtx<double>(data_dir + "time_mesh.mtx");
    int n_times = time_mesh.rows();
    int n_sim = 30;
    int n_locs = 250;
    std::vector<std::string> sigma_lab = {"0.00", "0.05", "0.10"};
    
    std::mt19937 gen(12345);
    for(int i=0; i<3; ++i){
    std::string sigma_dir = data_dir + "sigma_" + sigma_lab[i] + "/";
    for(int k = 0; k < n_sim; ++k){
        
        std::string sim_dir = sigma_dir + std::to_string(k) + "/";

        auto solution = read_mtx<double>(sim_dir + "exact.mtx");
    
        matrix_t locs = unit_square.sample(n_locs);

        auto psi = internals::point_basis_eval(Vh, locs);
        matrix_t obs = matrix_t::Zero(n_locs, n_times-1); // t0 NON prendo osservazioni (?)

        for(int t=0; t < n_times-1; ++t) obs.col(t) = psi * solution.col(t+1) + noise(n_locs, 0.05, gen);

        Eigen::saveMarket(locs, sim_dir + "locs.mtx");
        Eigen::saveMarket(obs, sim_dir + "obs.mtx");
        }
    }

    std::string command_str = "chown -R 1000:1000 " + data_dir; 
    system(command_str.c_str());
    
    return 0;
}
