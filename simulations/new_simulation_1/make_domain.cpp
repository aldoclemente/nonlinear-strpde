
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

    using FeType = FeP<1, 1>;
    using Quadrature = FeType::template cell_quadrature_t<local_dim>;

    Triangulation<2, 2> unit_square = Triangulation<2, 2>::UnitSquare(21, cache_cells);
    
    std::string mesh_dir = "data/mesh/";

    std::string command_str = "mkdir -p " + mesh_dir; 
    system(command_str.c_str());

    Eigen::saveMarket(unit_square.nodes(), mesh_dir + "points.mtx");
    Eigen::saveMarket(unit_square.cells(), mesh_dir + "cells.mtx");
    Eigen::saveMarket(unit_square.boundary_nodes().as_eigen_matrix(), mesh_dir + "boundary.mtx");
    
    auto quad_nodes = simplex_quadrature_nodes(unit_square, Quadrature{});
    Eigen::saveMarket(quad_nodes, mesh_dir + "quad_nodes.mtx");

    double T = 1.; 
    int n_intervals = 10;
    int n_times = n_intervals + 1;
    double DeltaT = T/n_intervals;
    
    vector_t time_mesh = vector_t::Zero(n_times);
    for(int t = 0; t < n_times; t++) time_mesh[t] = DeltaT*t;
    Eigen::saveMarket(time_mesh, mesh_dir + "time_mesh.mtx");

    command_str = "chown -R 1000:1000 " + mesh_dir; 
    system(command_str.c_str());

    return 0;

}   