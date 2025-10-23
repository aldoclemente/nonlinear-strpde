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

    std::string mesh_dir = "data/mesh/";

    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);
    
    int n_dofs = unit_square.nodes().rows();
    matrix_t test_locs = read_mtx<double>(mesh_dir + "test_locs.mtx");
    
    // FE space
    
    FeSpace Vh(unit_square, P1<1>);   // piecewise linear continuous scalar finite elements
    TrialFunction u(Vh);
    TestFunction  v(Vh);

    FeFunction u_prev(Vh);

    matrix_t IC = read_mtx<double>(mesh_dir + "IC.mtx");
    
    vector_t time_mesh = read_mtx<double>(mesh_dir + "time_mesh.mtx");
    time_mesh*= (2e-2);
    double DeltaT = time_mesh[1] - time_mesh[0];
    std::cout << DeltaT << std::endl;
    int n_times = time_mesh.size();
    std::cout << IC.rows() << std::endl;
    u_prev = IC;
    std::cout << "ciao" << std::endl;
    auto psi_test = internals::point_basis_eval(Vh, test_locs);
    
    // laplacian operator bilinear form
    double mu = 1.;
    int n_sim = 30;
    double r = 10.0;
    //matrix_t alpha_mtx = matrix_t::Ones(n_sim, 3);
    //alpha_mtx *= r;

    //std::vector<std::string> sigma_lab = {"0.00", "0.05", "0.10"};
    //std::mt19937 gen(12345);
    //for(int i=1; i < 3; ++i) alpha_mtx.col(i) += noise(n_sim, 0.05*i, gen);

    //std::string coeff_dir = data_dir + "mu_"  + std::to_string(mu) + "_r_" + std::to_string(r) + "/";
    //std::string command_str = "mkdir -p " + coeff_dir; 
    //system(command_str.c_str());
    
    auto a = integral(unit_square)(mu*dot(grad(u), grad(v)) - r * u * v);
    auto reac = integral(unit_square)(r * u_prev*u*v);
    auto mass = integral(unit_square)(u*v);

    auto f_coeff = read_mtx<double>(mesh_dir + "forcing_coeff.mtx");
    f_coeff *=50;
    matrix_t solution = matrix_t::Zero(n_dofs, n_times);
    matrix_t parabolic = matrix_t::Zero(n_dofs, n_times);

    solution.col(0) = IC;
    parabolic.col(0) = IC;
    sparse_matrix_t M = mass.assemble();
    sparse_matrix_t A = a.assemble();
    vector_t rhs = vector_t::Zero(n_dofs);

    //f_coeff ;
    FeCoeff<local_dim, 1, 1, vector_t> f(f_coeff.col(0));
    auto F = integral(unit_square)(f*v).assemble();
    std::cout << f_coeff.col(0).minCoeff() << " " << f_coeff.col(0).maxCoeff() << std::endl;
    std::cout << F.minCoeff() << " " << F.maxCoeff() << std::endl;
    for(int t = 0; t < n_times-1; ++t){
    // NONLINEAR ---------------------------------------------------------
    u_prev = solution.col(t);
    auto R = reac.assemble();
    
    auto S = (1./DeltaT*M + A + R);
    rhs = 1./DeltaT * M * solution.col(t) + F; // F(t+1)+

    sparsesolver_t lin_solver(S);
        
    lin_solver.factorize(S);
        
    solution.col(t+1) = lin_solver.solve(rhs);

    // PARABOLIC ---------------------------------------------------------

    auto S2 = (1./DeltaT*M + A);
    rhs = 1./DeltaT * M * parabolic.col(t) + F;
    lin_solver.factorize(S2);
    parabolic.col(t+1) = lin_solver.solve(rhs);
    }

    matrix_t test_obs = matrix_t::Zero(test_locs.rows(), n_times); // t0 butto via
    for(int t = 0; t<n_times; ++t) test_obs.col(t) = psi_test * solution.col(t);

    Eigen::saveMarket(solution, mesh_dir + "exact.mtx");

    Eigen::saveMarket(f_coeff.col(0), mesh_dir + "forcing_coeff_pde_sol.mtx");
    Eigen::saveMarket(time_mesh, mesh_dir + "time_mesh_pde_sol.mtx");
    Eigen::saveMarket(parabolic, mesh_dir + "parabolic.mtx");
    Eigen::saveMarket(test_obs, mesh_dir + "test_obs.mtx");

    std::string command_str = "chown -R 1000:1000 " + mesh_dir; 
    system(command_str.c_str());
    return 0;
}
    