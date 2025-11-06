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
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>::UnitSquare(31, cache_cells);
    
    std::string data_dir = "input/";
    std::string mesh_dir = "input/mesh/";
    
    std::string command_str = "mkdir -p " + mesh_dir; 
    system(command_str.c_str());
    
    Eigen::saveMarket(unit_square.nodes(), mesh_dir + "points.mtx");
    Eigen::saveMarket(unit_square.cells(), mesh_dir + "cells.mtx");
    Eigen::saveMarket(unit_square.boundary_nodes().as_eigen_matrix(), mesh_dir + "boundary.mtx");
    //Eigen::saveMarket(unit_square, mesh_dir + "quad_nodes.mtx");

    // FE space
    auto quad = P1<1>;
    
    FeSpace Vh(unit_square, quad);   // piecewise linear continuous scalar finite elements
    TrialFunction u(Vh);
    TestFunction  v(Vh);
    
    auto& dof_handler = Vh.dof_handler();
    int n_dofs = dof_handler.n_dofs();
    matrix_t dofs =  unit_square.nodes(); // fe_dofs_coords(unit_square, quad);
    std::cout << dofs.rows() << " " << dofs.cols() << std::endl;
    unit_square.mark_boundary(/* as = */ 0);
    unit_square.mark_boundary(/* as = */ 4, /* where = */ [](const typename Triangulation<2, 2>::EdgeType& e) {
       return (e.node(0)[0] == 0 && e.node(1)[0] == 0);
    });
    //dof_handler.set_dirichlet_constraint(/* on = */ 0, /* data = */ g_D);
    
    FeFunction u_prev(Vh);

    // define forcing functional
    
    // ScalarField<2, decltype([](const PointT& p) {
    // if        (  ((p[0]-0.5)*(p[0]-0.5) + (p[1]-0.5)*(p[1]-0.5)) <= 0.1*0.1 ) {
    //     return 1.;
    // } else {
    //     return 0.;
    // }})> ic;
    
    auto coe = [](double x) -> double {
        return 0.5*std::sin(5.*std::numbers::pi*x)*std::exp(-x*x) + 1.;
    };

    double tf = 1.0;
    double t = 0.0;
    
    auto vicini = [&tf, &t, &coe](const PointT& p) {
        double a = (9./(5.*tf))*t - 2.;
        return abs( std::sin( 2.*std::numbers::pi*(  coe(p[1])*p[0]*cos(a) - p[1]*std::sin(a) ) ) * 
                    std::cos( 2.*std::numbers::pi*(  coe(p[1])*p[0]*cos(a + std::numbers::pi/2.) + coe(p[0])*p[1]*std::sin( (std::numbers::pi/2.)*a ) ) )
        );
    };

    ScalarField<2> ic;
    ic = vicini;

    // ScalarField<2, decltype([](const PointT& p) {
    //     double sigma = 0.10;
    //     return std::exp( -((p[0]-0.5)*(p[0]-0.5) + (p[1]-0.5)*(p[1]-0.5)) / (2.*sigma*sigma) );
        
    // })> ic;
    
    
    double DeltaT = 1e-1;
    int n_times = 11;
    double T = (n_times-1) * DeltaT;
    vector_t time_mesh = vector_t::Zero(n_times);
    for(int t = 0; t < n_times; t++) time_mesh[t] = DeltaT*t;
    
    Eigen::saveMarket(time_mesh, mesh_dir + "time_mesh.mtx");

    vector_t IC = vector_t::Zero(n_dofs);
    for(int i = 0; i < n_dofs; ++i) IC[i] = ic(dofs.row(i));
    std::cout << dofs.rows() << " " << dofs.cols() << std::endl;
    std::cout << unit_square.nodes().rows() << " " << unit_square.nodes().cols() << std::endl;
    
   // vector_t IC = read_mtx<double>(mesh_dir + "IC.mtx");
    std::cout << IC.minCoeff() << " " << IC.maxCoeff() << std::endl;
        
    u_prev = IC;


    // laplacian operator bilinear form
    double mu = 1e-2; // 0.01
    double alpha = 1.;
    int n_sim = 30;

        //double alpha = alpha_mtx(k,i);
        auto a = integral(unit_square)(mu*dot(grad(u), grad(v)));
        auto reac= integral(unit_square)(u_prev*u*v);
        auto mass = integral(unit_square)(u*v);

        matrix_t solution = matrix_t::Zero(n_dofs, n_times);
        matrix_t diffusion = matrix_t::Zero(n_dofs, n_times);
        
        solution.col(0) = IC;
        diffusion.col(0) = IC;

        auto M = mass.assemble();
        auto A = a.assemble();
        vector_t rhs = vector_t::Zero(n_dofs);
        
        for(int t = 0; t < n_times-1; ++t){
            u_prev = solution.col(t);
            auto R = reac.assemble();
            range(R);
           
            //semi-implicit (no mass lumping) 
            sparse_matrix_t S = (1./DeltaT*M + A) - alpha*(M - R);
            std::cout << "S: " << std::endl; 
            range(S);
            
            rhs = 1./DeltaT * M*solution.col(t);
            
            sparsesolver_t lin_solver;
            lin_solver.analyzePattern(S); 
            lin_solver.factorize(S);
            solution.col(t+1) = lin_solver.solve(rhs);
            std::cout << "range f("<<t+1<<") " << solution.col(t + 1).minCoeff() << " " << solution.col(t + 1).maxCoeff() << "\n" << std::endl;

            //---
            S = (1./DeltaT*M + A);
            rhs = 1./DeltaT * M*diffusion.col(t);
            lin_solver.analyzePattern(S); 
            lin_solver.factorize(S);
            diffusion.col(t+1) = lin_solver.solve(rhs);  
        }

        matrix_t test_locs = expand_grid(std::vector<vector_t>{
                                         vector_t::LinSpaced(60, 0.05, 0.95), 
                                         vector_t::LinSpaced(60, 0.05, 0.95)});
        auto psi_test = internals::point_basis_eval(Vh, test_locs);
    
        matrix_t test_obs = matrix_t::Zero(test_locs.rows(), n_times); 
        for(int t = 0; t<n_times; ++t) test_obs.col(t) = psi_test * solution.col(t);

        Eigen::saveMarket(solution, data_dir + "fisher_kpp.mtx");
        Eigen::saveMarket(diffusion, data_dir + "diffusion.mtx");
        
        Eigen::saveMarket(test_obs, mesh_dir + "test_obs.mtx");
        Eigen::saveMarket(test_locs, mesh_dir + "test_locs.mtx");


    command_str = "chown -R 1000:1000 " + data_dir; 
    system(command_str.c_str());
    return 0;
}
