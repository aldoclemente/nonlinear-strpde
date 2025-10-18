
#include <LBFGSB.h>
#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;

Eigen::Matrix<double, Dynamic, 1> seq(double start, double end, double step) {
    // number of points (inclusive), robust to FP noise
    int n = static_cast<int>(std::floor((end - start) / step + 1e-12)) + 1;
    // clamp exact end to start + step*(n-1) to keep spacing exact
    double adjusted_end = start + step * (n - 1);
    return Eigen::Matrix<double, Dynamic, 1>::LinSpaced(n, start, adjusted_end);
}

int main(int argc, char *argv[]){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "../data/simulation_1/";
    std::string mesh_dir = "../data/simulation_1/mesh/";
    
    int sd = std::stoi(argv[1]); // 0 -> "0.00", 1 -> "0.05", 3 -> "0.10" 
    std::string sim = std::string(argv[2]); // 0, ..., 29

    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    matrix_t test_locs = read_mtx<double>(data_dir + "test_locs.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);

    matrix_t time_mesh = read_mtx<double>(data_dir + "time_mesh.mtx");
    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    int n_times = time_mesh.rows();
    int n_sim = 30;
    Triangulation<1, 1> T(time_mesh);

    std::vector<std::string> sigma_lab = {"0.00", "0.05", "0.10"};

    std::string sigma_dir = data_dir + "sigma_" + sigma_lab[sd] + "/";
    
    std::string out_dir = "../data/param_cascading_sim_1/sigma_" + sigma_lab[sd] + "/" + sim + "/";

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
    double alpha = 3.0; // true alpha
    
    vector_t rho_seq = vector_t::Zero(13, 1);
    rho_seq << 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99;

    //int n_rho = rho_seq.rows();
    //vector_t gcv_vec = vector_t::Zero(n_rho, 1);
    //matrix_t alpha_optim = matrix_t::Zero(n_rho, 2);

    // GCV (inner opt criterion)
    auto gamma = seq(-4.0,4.0,0.25);
    int n_lambda = gamma.rows();
    GridSearch<2> optimizer;
    matrix_t rho_grid = matrix_t::Ones(n_lambda, 2);
    for (int i = 0; i < n_lambda; ++i) { rho_grid(i, 0) = std::pow(10, gamma(i,0))/(1 + std::pow(10, gamma(i,0))); }
    
    //rho_grid.col(0) = vector_t::LinSpaced(n_lambda, 1e-3, 0.995);

    matrix_t lambda_grid = matrix_t::Ones(rho_grid.rows(), 2);
    lambda_grid.col(0) = rho_grid.col(0).array() / (unit_square.measure() * (1.0 - rho_grid.col(0).array()));
    
    std::cout << rho_grid.transpose() << std::endl;
    std::cout << lambda_grid.transpose() << std::endl;
    
    //for (int i = 0; i < n_lambda; ++i) { rho_grid(i, 0) = std::pow(10, -5.0 + 0.25 * i); }
    
    // initial guess
    vector_t x0 = vector_t::Ones(1);
    // SSE (outer opt criterion) 
    double rho = 1.;
    size_t i = 0;
    ScalarField<1> PC;
    PC = [&unit_square, &data, &IC, &mu, &lambda_grid, &i](vector_t x) -> double {
        FeSpace Vh(unit_square, P1<1>);
        TrialFunction u(Vh);
        TestFunction v(Vh);
        ZeroField<2> f;
        auto a = integral(unit_square)(mu * dot(grad(u), grad(v)) - x(0,0) * u * v);
        auto F = integral(unit_square)(f * v);
        SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));
        
        // inner fitting criterion (fixed lambda)
        parabolic.fit(lambda_grid.row(i));

        // outer opt criterion
        return ((parabolic.fitted() - parabolic.response()).squaredNorm());
    };
    
    struct function_to_optimize {
        function_to_optimize(const ScalarField<1>& Field) : Field_(Field) { }

        double operator()(const vector_t& x, vector_t& grad) {
            double Fx = Field_(x);
            grad = Field_.gradient()(x);
            return Fx;
        }

        const ScalarField<1>& Field_;
    };

    function_to_optimize F_(PC);

    // Set up parameters
    LBFGSpp::LBFGSBParam<double> param;
    param.max_iterations = 100;
    param.epsilon_rel = 1e-6;
    LBFGSpp::LBFGSBSolver<double> bfgsb(param);

    // theta, gamma
    vector_t lower_bounds = vector_t::Zero(1, 1);
    vector_t upper_bounds = vector_t::Zero(1, 1);
    upper_bounds << 10.;

    NelderMead<1> nelder_mead(30, 1e-4);
    //nelder_mead.optimize(PC, x0,BacktrackingLineSearch());
    vector_t alpha_optim = vector_t::Zero(n_lambda);
    vector_t gcv_vec = vector_t::Zero(n_lambda);

    auto GCV = [&unit_square, &data, &IC, &mu, &lambda_grid, &i, &alpha_optim, &optimizer] () -> double {
        FeSpace Vh(unit_square, P1<1>);
        TrialFunction u(Vh);
        TestFunction v(Vh);
        ZeroField<2> f;
        
        auto a = integral(unit_square)(mu * dot(grad(u), grad(v)) - alpha_optim(i, 0) * u * v);
        auto F = integral(unit_square)(f * v);
        SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));

        //optimizer.optimize(parabolic.gcv(100, 476813), lambda_grid);
        
        //parabolic.gcv(100, 476813)(lambda_grid.row(i))
        
        //parabolic.fit(optimizer.optimum());
        return parabolic.gcv(100, 476813)(lambda_grid.row(i));
    };

    // for(size_t i=0; i< rho_grid.rows(); ++i){
    //     std::cout << "\t --- rho " << rho_grid(i, 0) << " ---" << std::endl;
    //     vector_t x = x0;
    //     if (i != 0) x(0, 0) = alpha_optim(i - 1, 0);
    
    //     nelder_mead.optimize(PC,x, BacktrackingLineSearch());
    //     std::cout << "NELDER MEAD" << std::endl;
    //     std::cout << "x opt: " << nelder_mead.optimum()[0]  << std::endl;
    //     std::cout << "value opt: " << nelder_mead.value() << std::endl;
    //     std::cout << "#iters: " << nelder_mead.n_iter() << std::endl;
    //     alpha_optim(i, 0) = nelder_mead.optimum()[0];

    //     gcv_vec(i, 0) = GCV();
    // }
    
    for(size_t i=0; i< rho_grid.rows(); ++i){
        std::cout << "\t --- rho " << rho_grid(i, 0) << " ---" << std::endl;
        vector_t x = x0;
        if (i != 0) x(0, 0) = alpha_optim(i - 1, 0);

        double fx;
        int niter = bfgsb.minimize(F_, x, fx, lower_bounds, upper_bounds);
        std::cout << "\t LBFGS" << std::endl;
        std::cout << "iter: " << niter << std::endl;
        std::cout << "optimum: " << x.transpose() << std::endl;
        alpha_optim(i, 0) = nelder_mead.optimum()[0];

        gcv_vec(i, 0) = GCV();
    }

    std::cout << "\ngcv_vec: \n" << gcv_vec << std::endl;
    // Best alpha value?!
    Eigen::Index minRow;
    gcv_vec.col(0).minCoeff(&minRow);
    std::cout << "gcv min at " << minRow << "  value: " << gcv_vec(minRow, 0) << std::endl;
    
    std::string command_str = "mkdir -p " + out_dir;
    system(command_str.c_str());

    //Eigen::saveMarket(parabolic.f(), sim_dir + "estimate_parabolic.mtx");

    // rmse
    //auto psi_test = internals::point_basis_eval(Vh, test_locs);

    //int n_dofs = unit_square.nodes().rows();
    //matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
    // for(int t = 0; t<n_times-1; ++t){ 
    //     test_vals.col(t) = psi_test * parabolic.f().block(n_dofs*t, 0, n_dofs, 1);    
    // }

    // matrix_t test_obs = read_mtx<double>(sim_dir + "test_obs.mtx");
    
    // rmse[0] = std::sqrt( (test_vals - test_obs).reshaped().array().square().mean()); 
    
    // Eigen::saveMarket(rmse, sim_dir + "rmse_parabolic.mtx");
    
    // command_str = "chown -R 1000:1000 " + data_dir; 
    // system(command_str.c_str());

    return 0;
}