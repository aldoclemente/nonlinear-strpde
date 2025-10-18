#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
#include "../include/fe_ls_fisher_kpp.h"
#include "../include/kFoldCV.h"
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;

int main(int argc, char *argv[]){
    std::cout << "\t sim 3 init" << std::endl;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "../data/simulation_3/";
    std::string data_dir2 = "../data/simulation_2/";
    std::string mesh_dir = "../data/simulation_1/mesh/";
    
    int sd = std::stoi(argv[1]); // 0 -> "0.00", 1 -> "0.05", 3 -> "0.10" 
    std::string sim = std::string(argv[2]); // 0, ..., 29

    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    
    matrix_t test_locs = read_mtx<double>("../data/simulation_1/test_locs.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);
    std::vector<std::string> sigma_lab = {"0.00", "0.05", "0.10"};
    
    FeSpace Vh(unit_square, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);
    ZeroField<local_dim> f;

    Eigen::Matrix<int, Dynamic, Dynamic> incidence_matrix = read_csv<int>(data_dir + "incidence_matrix.csv").as_matrix();
    

    std::string sigma_dir = data_dir + "sigma_" + sigma_lab[sd] + "/";
    std::string sigma_dir2 = data_dir2 + "sigma_" + sigma_lab[sd] + "/";

    std::string sim_dir = sigma_dir + sim + "/";
    std::string sim_dir2 = sigma_dir2 + sim + "/";
    
    vector_t obs = read_mtx<double>(sim_dir + "obs.mtx").col(0);
    
    double mu = 1.0;

    GeoFrame data(unit_square);
    auto bm = BinaryMatrix<Dynamic, Dynamic> (incidence_matrix.rows(), incidence_matrix.cols());
    bm = incidence_matrix;

    std::cout << "incidence " << incidence_matrix.rows() << " " << incidence_matrix.cols() << std::endl;
    std::cout << "BM " << bm.rows() << " " << bm.cols() << std::endl;
    std::cout << ":=)" << std::endl;
    auto &l = data.insert_scalar_layer<POLYGON>(
                "layer", bm); 
    l.load_vec("y", obs.reshaped());

    auto a = integral(unit_square)(mu * dot(grad(u), grad(v)));
    auto F = integral(unit_square)(f * v);
    SRPDE model("y ~ f", data, fe_ls_elliptic(a, F));

    int n_lambda = 50;
    matrix_t lambda_grid = matrix_t::Ones(n_lambda+1,1);
    for(int i=0; i<=n_lambda;++i) lambda_grid(i,0) = std::pow(10, -5.0 + 0.1 * i);

    GridSearch<1> optimizer;
    /* 
    std::cout << ": = )" << std::endl;
    optimizer.optimize(model.gcv(100, 476813), lambda_grid);
    // ???????? rotto ??????
    std::cout << ":=( )" << std::endl;

    std::cout << "lambda " << optimizer.optimum() << std::endl;
    Eigen::saveMarket(optimizer.optimum(), data_dir + "IC_lambda.mtx");

    std::cout << "fit " << std::endl;
    
    model.fit(optimizer.optimum()[0]);

    matrix_t values = matrix_t::Zero(n_lambda+1,1);
    for(int i=0; i<=n_lambda; ++i) values(i,0) = optimizer.values()[i];
    Eigen::saveMarket(values, data_dir + "IC_gcv.mtx"); */

    double lambda = 1e-3;
    std::cout << "- fit" << std::endl;
    model.fit(lambda);
    std::cout << model.f().transpose() << std::endl;
    //std::cout << "- save output" << std::endl;
    //Eigen::saveMarket(model.f(), data_dir + "IC.mtx");

    //std::string command_str = "chown -R 1000:1000 " + data_dir; 
    //system(command_str.c_str());
    std::cout << "\t ended" << std::endl;
    return 0;
}