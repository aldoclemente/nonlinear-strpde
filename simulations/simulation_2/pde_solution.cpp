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
    
    int n_dofs = unit_square.nodes().rows();
    matrix_t test_locs = read_mtx<double>("../data/simulation_1/test_locs.mtx");
    
    // FE space
    
    FeSpace Vh(unit_square, P1<1>);   // piecewise linear continuous scalar finite elements
    TrialFunction u(Vh);
    TestFunction  v(Vh);

    FeFunction u_prev(Vh);

    matrix_t IC = read_mtx<double>(data_dir + "ic.mtx");
    
    double T = 0.5, DeltaT = 0.05;
    int n_times = std::ceil(T/DeltaT) + 1; 
    vector_t time_mesh = vector_t::Zero(n_times);
    for(int t = 0; t < n_times; t++) time_mesh[t] = DeltaT*t;
    Eigen::saveMarket(time_mesh, data_dir + "time_mesh.mtx");

    u_prev = IC;

    auto psi_test = internals::point_basis_eval(Vh, test_locs);
    
    // laplacian operator bilinear form
    double mu = 0.01;
    int n_sim = 30;

    matrix_t alpha_mtx = matrix_t::Ones(n_sim, 3);
    alpha_mtx *= 3.0;

    std::vector<std::string> sigma_lab = {"0.00", "0.05", "0.10"};
    std::mt19937 gen(12345);
    for(int i=1; i < 3; ++i) alpha_mtx.col(i) += noise(n_sim, 0.05*i, gen);

    for(int i=0; i<3; ++i){
    std::string sigma_dir = data_dir + "sigma_" + sigma_lab[i] + "/";
    
    if (!std::filesystem::exists(
        std::filesystem::path(sigma_dir))) {
            std::filesystem::create_directory(
                std::filesystem::path(sigma_dir));
    }

    for(int k=0; k<n_sim; ++k){

        std::string sim_dir = sigma_dir + std::to_string(k) + "/";
        if (!std::filesystem::exists(
            std::filesystem::path(sim_dir))) {
            std::filesystem::create_directory(
                std::filesystem::path(sim_dir));
        }

        double alpha = alpha_mtx(k,i);
        auto a = integral(unit_square)(mu*dot(grad(u), grad(v)) - alpha * u * v);
        auto reac = integral(unit_square)(alpha * u_prev*u*v);
        auto mass = integral(unit_square)(u*v);

        matrix_t solution = matrix_t::Zero(n_dofs, n_times);
        solution.col(0) = IC;
        auto M = mass.assemble();
        auto A = a.assemble();
        vector_t rhs = vector_t::Zero(n_dofs);

        for(int t = 0; t < n_times-1; ++t){

            u_prev = solution.col(t);
            auto R = reac.assemble();

            auto S = (1./DeltaT*M + A + R);
            rhs = 1./DeltaT * M * solution.col(t);

            sparsesolver_t lin_solver(S);
        
            lin_solver.factorize(S);
        
            solution.col(t+1) = lin_solver.solve(rhs);
        }

        matrix_t test_obs = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
        for(int t = 0; t<n_times-1; ++t) test_obs.col(t) = psi_test * solution.col(t+1);

        Eigen::saveMarket(solution, sim_dir + "exact.mtx");
        Eigen::saveMarket(test_obs, sim_dir + "test_obs.mtx");
        }
    }

    std::string command_str = "chown -R 1000:1000 " + data_dir; 
    system(command_str.c_str());
    return 0;
}
