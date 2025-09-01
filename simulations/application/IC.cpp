#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
#include "../include/fe_ls_fisher_kpp.h"
#include "../include/kFoldCV.h"
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;

int main(){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 3;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "../data/application/";
    std::string mesh_dir = data_dir + "mesh/";

    matrix_t nodes = read_csv<double>(mesh_dir + "nodes.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> cells = read_csv<int>(mesh_dir + "cells.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>(mesh_dir + "boundary.csv").as_matrix();

    Triangulation<3, 3> brain = Triangulation<3, 3>(nodes, cells, boundary);

    FeSpace Vh(brain, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);
    ZeroField<3> f;

    Eigen::Matrix<int, Dynamic, Dynamic> incidence_matrix = read_csv<int>(mesh_dir + "incidence_matrix.csv").as_matrix();
    
    matrix_t obs = read_csv<double>(data_dir + "obs0.csv").as_matrix();
    
    double mu = 1.0;

    GeoFrame data(brain);
    auto bm = BinaryMatrix<Dynamic, Dynamic> (incidence_matrix.rows(), incidence_matrix.cols());
    auto &l = data.insert_scalar_layer<POLYGON>(
                "layer", bm); 
    l.load_vec("y", obs.reshaped());

    auto a = integral(brain)(mu * dot(grad(u), grad(v)));
    auto F = integral(brain)(f * v);
    SRPDE model("y ~ f", data, fe_ls_elliptic(a, F));

    int n_lambda = 50;
    matrix_t lambda_grid = matrix_t::Ones(n_lambda+1,1);
    for(int i=0; i<=n_lambda;++i) lambda_grid(i,0) = std::pow(10, -5.0 + 0.1 * i);

    GridSearch<1> optimizer;
    optimizer.optimize(model.gcv(100, 476813), lambda_grid);
    
    model.fit(optimizer.optimum());
    Eigen::saveMarket(model.f(), data_dir + "IC.mtx");
}