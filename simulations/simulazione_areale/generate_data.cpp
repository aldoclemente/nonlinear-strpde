#include <fdaPDE/core/fdaPDE/core.h>
#include <fdaPDE/src/solvers/utility.h>
#include "../include/utils.h"
#include <unsupported/Eigen/SparseExtra>
#include <string>
#include <vector>

using namespace fdapde;
int main (){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "input/";
    std::string mesh_dir = "input/mesh/";
    
    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);
    FeSpace Vh(unit_square, P1<1>);   // piecewise linear continuous scalar finite elements
    
    std::string command_str;
    vector_t time_mesh = read_mtx<double>(mesh_dir + "time_mesh.mtx");
    int n_times = time_mesh.size();
    std::cout << n_times << std::endl;
    std::cout << time_mesh[1] - time_mesh[1] << std::endl; 
    std::cout << time_mesh << std::endl;
    int n_sim = 30;
    
    auto solution = read_mtx<double>(data_dir + "fisher_kpp.mtx");
    std::mt19937 gen(12345);

    Eigen::Matrix<int, Dynamic, Dynamic> incidence_matrix = read_csv<int>(data_dir + "incidence_matrix.csv").as_matrix();
    auto bm = BinaryMatrix<Dynamic, Dynamic> (incidence_matrix.rows(), incidence_matrix.cols());
    bm = incidence_matrix;

    const auto &[psi, measure_vec] =
        internals::areal_basis_eval(Vh, bm);
    
    for(int i=0; i < n_sim; ++i){
        std::string sim_dir = data_dir + std::to_string(i) + "/";
        
        matrix_t obs = matrix_t::Zero(bm.rows(), n_times); 
            
        for(int t=0; t < n_times; ++t) obs.col(t) = psi * solution.col(t);
        
        Eigen::saveMarket(obs, sim_dir + "obs_no_noise.mtx");
        obs += noise(bm.rows(), n_times, 0.05, gen);

        command_str = "mkdir -p " + sim_dir; 
        system(command_str.c_str());
    
        Eigen::saveMarket(obs, sim_dir + "obs.cpp.mtx");
    }

    return 0;
}