
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
    
    matrix_t obs = read_csv<double>(data_dir + "obs.csv").as_matrix();
    
    double mu = 1.0;

    matrix_t time_mesh = matrix_t::Zero(3,1);
    time_mesh(1,0) = 1.0;
    time_mesh(2,0) = 2.0;

    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);

    int n_times = time_mesh.rows();
    Triangulation<1, 1> T(time_mesh);

    // 10-fold cv to select the best model
    int n_folds = 10; 
    KFoldCV cv_outer(incidence_matrix, obs, n_folds);

    // first col: nonlin
    // second col: parabolic (without alpha)
    // third col : parabolic (with alpha)
    matrix_t cv_error = matrix_t::Zero(n_folds,3);

    double alpha = -0.143; // see Shafer et al.
    vector_t IC = read_mtx<double>(data_dir + "IC.mtx");

    for(int iter=0; iter<10; ++iter){
        
        auto incidence_matrix_iter = cv_outer.X_train(iter);
        auto obs_iter = cv_outer.Y_train(iter);
    
        GeoFrame data(brain, T);
        auto bm = BinaryMatrix<Dynamic, Dynamic> (incidence_matrix_iter.rows(), incidence_matrix_iter.cols());
        auto &l = data.insert_scalar_layer<POLYGON, POINT>(
                "layer", std::pair{bm, time_locs}); 
        l.load_vec("y", obs_iter.reshaped()); 

        SRPDE nonlinear("y~f", data, fe_ls_fisher_kpp(std::pair(mu,alpha), Vh, IC)); 
        LBFGS<Dynamic> lbfgs;
        int n_dofs = nodes.rows();
        vector_t g_init = 0.5*(vector_t::Random(n_dofs*time_locs.size()) + vector_t::Ones(n_dofs*time_locs.size()));
        //nonlinear.fit(1.0, g_init, lbfgs, BacktrackingLinesearch());
    
        // 10-fold CV to select lambda ----------------------------------------------------------------------------------------
        int K = 10; 
        KFoldCV cv_(incidence_matrix_iter, obs_iter, K);
        vector_t cv_error = vector_t::Ones(K);

        ScalarField<1> SSE;

        SSE = [&brain, &T, &cv_, &time_locs, &IC, &Vh, &mu, &alpha, &n_times, &n_dofs, &g_init](vector_t lambda) -> double {
            if (lambda[0] <= 0.0) return std::numeric_limits<double>::max();
            double res = 0;
            LBFGS<Dynamic> lbfgs;
            for(int k=0; k < cv_.k(); ++k){
            
                auto incidence_mtx = cv_.X_train(k);
                auto obs = cv_.Y_train(k);

                GeoFrame data(brain, T);
                auto bm = BinaryMatrix<Dynamic, Dynamic> (incidence_mtx.rows(), incidence_mtx.cols());
                auto &l = data.insert_scalar_layer<POLYGON, POINT>(
                                "layer", std::pair{bm, time_locs}); 
                l.load_vec("y", obs.reshaped() ); 
            
                SRPDE model("y~f", data, fe_ls_fisher_kpp(std::pair(mu,alpha), Vh, IC)); 
            
                //vector_t g_init = 0.5*(vector_t::Random(n_dofs*time_locs.size()) + vector_t::Ones(n_dofs*time_locs.size()));
                model.fit(lambda[0], g_init, lbfgs); // BacktrackingLineSearch()

                Eigen::Matrix<int, Dynamic, Dynamic> test_locs = cv_.X_test(k);
                matrix_t test_obs = cv_.Y_test(k);
                auto test_bm = BinaryMatrix<Dynamic, Dynamic> (test_locs.rows(), test_locs.cols()); // 1 x ...
                test_bm = test_locs; 
                auto [psi_test, measure_vect] = internals::areal_basis_eval(Vh, test_bm);

                matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
                for(int t = 0; t<n_times-1; ++t){ 
                    test_vals.col(t) = psi_test * model.f().block(n_dofs*t, 0, n_dofs, 1);    
                }

                res += std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());
            }
            return 1./cv_.k() * res;
        };

        int n_lambda = 11;
        vector_t lambda_grid = vector_t::Ones(n_lambda);
        for(int i=0; i<n_lambda;++i) lambda_grid(i,0) = std::pow(10, -1. + 0.1 * i);
    
        GridSearch<1> optimizer;
        
        optimizer.optimize(SSE, lambda_grid); 
        nonlinear.fit(optimizer.optimum()[0], g_init, lbfgs);
        
        Eigen::Matrix<int, Dynamic, Dynamic> test_locs = cv_outer.X_test(iter);
        matrix_t test_obs = cv_outer.Y_test(iter);
        auto test_bm = BinaryMatrix<Dynamic, Dynamic> (test_locs.rows(), test_locs.cols()); // 1 x ...
        test_bm = test_locs; 

        auto [psi_test, measure_vect] = internals::areal_basis_eval(Vh, test_bm);
        matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
        for(int t = 0; t<n_times-1; ++t){ 
            test_vals.col(t) = psi_test * nonlinear.f().block(n_dofs*t, 0, n_dofs, 1);    
        }

        cv_error(iter,0) = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());
        // ----------------------------------------------------------------------------------------------------------------
        // PARABOLIC 1
        auto a_par = integral(brain)(mu * dot(grad(u), grad(v)));
        auto F_par = integral(brain)(f * v);
        SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a_par, F_par}, IC));

        optimizer.optimize(parabolic.gcv(100, 476813), lambda_grid);
        parabolic.fit(optimizer.optimum());

        test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
        for(int t = 0; t<n_times-1; ++t){ 
            test_vals.col(t) = psi_test * parabolic.f().block(n_dofs*t, 0, n_dofs, 1);    
        }

        cv_error(iter,1) = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());

        // PARABOLIC 2
        auto a_par2 = integral(brain)(mu * dot(grad(u), grad(v)) - alpha*u*v);
        SRPDE parabolic2("y ~ f", data, fe_ls_parabolic_mono(std::pair{a_par2, F_par}, IC));

        optimizer.optimize(parabolic2.gcv(100, 476813), lambda_grid);
        parabolic2.fit(optimizer.optimum());

        test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
        for(int t = 0; t<n_times-1; ++t){ 
            test_vals.col(t) = psi_test * parabolic2.f().block(n_dofs*t, 0, n_dofs, 1);    
        }

        cv_error(iter,2) = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean());

        Eigen::saveMarket(nonlinear.f(), data_dir  + "estimate_nonlinear_" +  std::to_string(iter) + ".mtx");
        Eigen::saveMarket(parabolic.f(), data_dir  + "estimate_parabolic_" +  std::to_string(iter) + ".mtx");
        Eigen::saveMarket(parabolic2.f(), data_dir + "estimate_parabolic2_"+  std::to_string(iter) + ".mtx");
    }
    
    Eigen::saveMarket(cv_error, data_dir  + "cv_error.mtx");
    std::string command_str = "chown -R 1000:1000 " + data_dir; 
    system(command_str.c_str());
    return 0;
}