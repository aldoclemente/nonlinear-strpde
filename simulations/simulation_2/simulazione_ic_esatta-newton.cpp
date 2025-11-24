#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
//#include "../include/fe_ls_fisher_kpp_OLD.h"
//#include "../include/fe_ls_fisher_kpp_lumped.h"
#include "../include/kFoldCV.h"
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;

struct fe_fisher_kpp {

private:
  void init() {

    FeSpace Vh(domain_, P1<1>);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);

    auto a =
        integral(domain_)(mu_ * dot(grad(uh), grad(vh)) - alpha_ * uh * vh);
    auto mass = integral(domain_)(uh * vh);

    A_ = a.assemble();
    A_.makeCompressed();

    M_ = mass.assemble();
    M_.makeCompressed();

    Psi_ = internals::point_basis_eval(a.trial_space(), locations_);
    Psi_.makeCompressed();

    PsiTPsi_ = Psi_.transpose() * Psi_;
    PsiTPsi_.makeCompressed();

    n_dofs_ = Vh.n_dofs();

    y_ = VectorType::Zero(n_dofs_ * n_time_locs_);
    u_ = y_;
    p_ = y_;

    SparseMatrixType I_(n_time_locs_, n_time_locs_);
    I_.setIdentity();

    Mt_ = kronecker(I_, M_);
    Mt_.makeCompressed();

    Psi_t_ = kronecker(I_, Psi_);

  }
  
public:
  using VectorType = Eigen::Matrix<double, Dynamic, 1>;
  using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
  using SparseMatrixType = Eigen::SparseMatrix<double>;
  using DiagonalMatrixType = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
  using SparseSolverType = Eigen::SparseLU<SparseMatrixType>;
  using DenseSolverType = Eigen::PartialPivLU<MatrixType>;

  fe_fisher_kpp(const Triangulation<2, 2> &domain, double mu, double alpha,
                const MatrixType &locations, const VectorType &time_locations,
                const VectorType &observations)
      : domain_(domain), mu_(mu), alpha_(alpha), locations_(locations),
        time_locations_(time_locations),
        dT_(time_locations[1] - time_locations[0]), observations_(observations),
        n_locs_(locations.rows()), n_time_locs_(time_locations.size()) {
    init();
  }

  void solve(double lambda) {
    lambda_ = lambda;
    u_ = opt_.optimize(*this, u_init_, BacktrackingLineSearch());
    n_iter_ = opt_.n_iter();
    y_ = state(u_);
    p_ = adjoint(y_);
  }

  double J(const VectorType &y, const VectorType &u) const {

    double SSE = 0;
    for (int t = 0; t < n_time_locs_; ++t) {
      SSE += (observations_.block(n_locs_ * t, 0, n_locs_, 1) - 
                Psi_ * y.block(n_dofs_ * t, 0, n_dofs_, 1)).squaredNorm(); //get_t(y, t)
    }
    return 1. / (n_locs_ * n_time_locs_) * SSE + lambda_ * u.squaredNorm();
  }

  VectorType state(const VectorType &u) const {

    FeSpace Vh(domain_, P1<1>);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    VectorType y = VectorType::Zero(n_dofs_ * n_time_locs_);
    y.block(0, 0, n_dofs_, 1) = y0_;

    for (int t = 0; t < n_time_locs_ - 1; t++) {
      y_old = y.block(n_dofs_ * t, 0, n_dofs_, 1);//get_t(y, t);
      auto R_ = reac.assemble();
      R_.makeCompressed();

      SparseMatrixType S = 1. / dT_ * M_ + A_ + R_;
      S.prune(1e-10);
      S.makeCompressed();

      SparseSolverType lin_solver(S);
      lin_solver.factorize(S);

      VectorType b = 1. / dT_ * M_ * y_old.coeff() + M_ *u.block(n_dofs_*(t+1), 0, n_dofs_, 1); //get_t(u, t + 1);

      y.block((t + 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return y;
  }

  VectorType adjoint(const VectorType &y) const {

    FeSpace Vh(domain_, P1<1>);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    VectorType p = VectorType::Zero(n_dofs_ * n_time_locs_); // p(T) == 0

    for (int t = n_time_locs_ - 1; t > 0; t--) {
      y_old = y.block(n_dofs_ * t, 0, n_dofs_, 1); 
      auto R_ = reac.assemble();
      R_.makeCompressed();

      SparseMatrixType S = 1. / dT_ * M_ + A_ + 2. * R_;
      S.prune(1e-10);
      S.makeCompressed();

      SparseSolverType lin_solver(S);
      lin_solver.factorize(S);

      VectorType b =
          1. / dT_ * M_ * p.block(n_dofs_ * t, 0, n_dofs_, 1) -
          1. / (n_time_locs_ * n_locs_) * PsiTPsi_ *y.block(n_dofs_*(t-1), 0, n_dofs_, 1) +
          1. / (n_time_locs_ * n_locs_) * Psi_.transpose() * observations_.block(n_locs_ * (t-1), 0, n_locs_, 1);

      p.block((t - 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return p;
  }

  double operator()(const VectorType &u) const {

    auto y = state(u);
    return J(y, u);
  }

  std::function<VectorType(const VectorType&)> gradient() {
    return [this](const VectorType &u) {
      auto y = state(u);
      auto p = adjoint(y);
      auto grad_u = lambda_ * Mt_ * u - Mt_ * p; // qua ci va lambda e non alpha (ci andrebbe un dT...)
      return grad_u;
    };
  }

  VectorType y() const { return y_; }

  VectorType u() const { return u_; }

  VectorType p() const { return p_; }

  VectorType get_t(const VectorType &v, int t) const {
    return v.block(n_dofs_ * t, 0, n_dofs_, 1);
  }

  //VectorType y(int t) const { return get_t(y(), t); }
  //VectorType u(int t) const { return get_t(u(), t); }
  //VectorType p(int t) const { return get_t(p(), t); }
  //VectorType obs(int t) const {
//     return observations_.block(n_locs_ * t, 0, n_locs_, 1);
//   }

  double J() const { 
    
    double res = (observations_ - Psi_t_ * y_).squaredNorm() + lambda_ * u_.squaredNorm();
    return res;  
  }

  int n_iter() const { return n_iter_;}

  // setters
  void set_state_initial_condition(const VectorType &y0) { y0_ = y0; }

  void set_control_initial_guess(const VectorType &u0) { u_init_ = u0; }

  // optimization algorithm custom stopping criterion
  template <typename OptimizerType> bool stop_if(OptimizerType &opt) {

    bool stop;

    VectorType y_old = state(opt.x_old);
    double loss_old = J(opt.x_old, y_old);

    VectorType y_new = state(opt.x_new);
    double loss_new = J(opt.x_new, y_new);

    stop = std::abs((loss_new - loss_old) / loss_old) < tol_;
    return stop;
  }

  // attribute
  const Triangulation<2, 2> &domain_;

  double mu_, alpha_;
  SparseMatrixType A_; // [A]_{ij} = \int_D mu * grad(Psi_j)*grad(Psi_i) - alpha
                       // * Psi_j * Psi_i
  SparseMatrixType M_; // Mass

  SparseMatrixType Psi_;
  SparseMatrixType PsiTPsi_;
  SparseMatrixType Mt_; // Spatio-Temporal Mass matrix
  SparseMatrixType Psi_t_;

  int n_dofs_ = 0;
  int n_locs_ = 0;
  int n_time_locs_ = 0;
  double dT_ = 0.;

  const MatrixType &locations_;
  const VectorType &time_locations_;
  VectorType observations_;

  double lambda_ = 1.;

  VectorType y_;
  VectorType u_;
  VectorType p_;

  VectorType y0_;
  VectorType u_init_;
  double tol_ = 1e-5;
  int n_iter_;
  LBFGS<Dynamic> opt_;
    
};

struct newton_method{

  using vector_t = Eigen::Matrix<double, Dynamic, 1>;
  using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
  using sparse_matrix_t = Eigen::SparseMatrix<double>;
  using sparse_solver_t = Eigen::SparseLU<sparse_matrix_t>;
  using dense_solver_t = Eigen::PartialPivLU<matrix_t>;
  
  void init() {

    FeSpace Vh(domain_, P1<1>);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);

    auto a =
        integral(domain_)(mu_ * dot(grad(uh), grad(vh)) - r_ * uh * vh);
    auto mass = integral(domain_)(uh * vh);

    A_ = a.assemble();
    A_.makeCompressed();

    M_ = mass.assemble();
    M_.makeCompressed();

    psi = internals::point_basis_eval(a.trial_space(), locations_);
    psi.makeCompressed();

    n_dofs = Vh.n_dofs();

    y_ = vector_t::Zero(n_dofs * m);
    u_ = y_;
    p_ = y_;

    sparse_matrix_t Im(m, m);
    Im.setIdentity();

    Mt_ = kronecker(Im, M_);
    Mt_.makeCompressed();

    Psi = kronecker(Im, psi);

    At_ = kronecker(Im, A_);
   
    // assemble matrix associated with derivation in time L_
    // [L]_{ii} = 1/DeltaT for i \in {1 ... m} and [L_]_{i,i-1} = -1/DeltaT for i \in {1 ... m-1}
    std::vector<Eigen::Triplet<double>> triplet_list;
    triplet_list.reserve(2 * m);
    // start assembly loop
    triplet_list.emplace_back(0, 0, 1.0 / deltaT );
    for (int i = 1; i < m; ++i) {
        triplet_list.emplace_back(i, i,  1.0 / deltaT );
        triplet_list.emplace_back(i, i - 1, -1.0 / deltaT );
    }
    sparse_matrix_t L(m, m);   // m x m identity matrix
    L.setFromTriplets(triplet_list.begin(), triplet_list.end());
    L.makeCompressed();
    D = kronecker(L, M_);

    Mw.resize(m);
    Gw.resize(m);
  }


  void update_weighted_mass(const vector_t& f_k, std::vector<sparse_matrix_t>& Mw) {
    FeSpace Vh(domain_, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);

    auto f_old = FeFunction(Vh, vector_t::Zero(n_dofs));
    auto mass_weight = integral(domain_)(f_old * u * v);
    
    for(int t=0; t<m; ++t){
      f_old = f_k.block(n_dofs * t, 0, n_dofs, 1);
      Mw[t] = mass_weight.assemble();
      Mw[t].makeCompressed();
      }
    };


  newton_method(const Triangulation<2, 2> &domain, double mu, double r,
                const matrix_t &locations, const vector_t &time_locations,
                const vector_t &observations)
      : domain_(domain), mu_(mu), r_(r), locations_(locations),
        time_locations_(time_locations),
        deltaT(time_locations[1] - time_locations[0]), observations_(observations),
        n(locations.rows()), m(time_locations.size()) {
      init();
        }
  
  void set_state_initial_condition(const vector_t &y0) { IC = y0; }

  std::pair<vector_t, int> solve(double lambda) {
        
        lambda_ = lambda;
        vector_t f  = vector_t::Zero(n_dofs*m);
        vector_t g  = vector_t::Zero(n_dofs*m);
        vector_t df = vector_t::Zero(n_dofs*m);
        vector_t dg = vector_t::Zero(n_dofs*m);
        // Sx = b
        sparse_solver_t invS_;
        vector_t dx = vector_t::Zero(2 * n_dofs * m);
        vector_t x = vector_t::Zero(2 * n_dofs * m);
        x.head(m*n_dofs) = f0;
        x.tail(m*n_dofs) = g0;

        vector_t b = vector_t::Zero(2 * n_dofs * m); // rhs 
        
        bool stop = false;
        int iter = 0;
        double tol = 1e-4;
        
        auto PsiTPsi = Psi.transpose()*Psi;
        
        update_weighted_mass(f0, Mw); 
        update_weighted_mass(g0, Gw); 
        
        vector_t u_ = vector_t::Zero(n_dofs * m); //no forcing
        u_.segment(0, n_dofs) += (1.0 / deltaT) * (M_*IC);
        double loss_old = std::numeric_limits<double>::max();
        double loss_new = 0;
        loss_old = J(x);
        //std::cout << "J: " << loss_old << std::endl;
        
        double sse_old = std::numeric_limits<double>::max();
        double sse_new = 0; 
        //std::cout << observations_.rows() << " " << observations_.cols() << std::endl;
        while( !stop && iter < max_iter){
            std::cout << "\t iter " << iter << std::endl;
            update_weighted_mass(f0, Mw); 
            update_weighted_mass(g0, Gw); 

            sparse_matrix_t A_tilde = blockDiag(Mw);
            //std::cout << A_tilde.rows() << " " << A_tilde.cols() << std::endl;
            //range(A_tilde);

            sparse_matrix_t G_tilde = blockDiag(Gw);
            //range(G_tilde);
            
            //std::cout << G_tilde.rows() << " " << G_tilde.cols() << std::endl;
            //std::cout << Psi.rows() << " " << Psi.cols() << std::endl;

            sparse_matrix_t S_00 = 1./(n * m)*Psi.transpose()*Psi + 2*r_*lambda*G_tilde;
            sparse_matrix_t S_10 = (D  + At_ -r_*Mt_ + 2*r_*A_tilde);

            sparse_matrix_t S_01 = lambda*S_10.transpose();

            SparseBlockMatrix<double, 2, 2> S ( S_00, S_01,
                                                S_10, -Mt_);
            
            invS_.compute(S);
            
            // right hand side
            b = vector_t::Zero(2 * n_dofs * m);

            b.head(n_dofs * m) = 1./(n * m) * Psi.transpose()*observations_ - 1./(n * m)* Psi.transpose()*(Psi*f0) - lambda*S_01*g0 ;
            b.tail(n_dofs * m) = -S_10*f0 + Mt_*g0;
            b.tail(n_dofs * m) += u_; 

            //std::cout << u_.rows() << " " << u_.cols() << std::endl;
            // solve
            dx = invS_.solve(b);
            df = dx.head(n_dofs * m);
            dg = dx.tail(n_dofs * m);

            //std::cout << "J_new: " << J(x+dx) << std::endl;
            //double alpha_opt = backtracking(x, dx);
            
            // std::cout << "J (backtracking): " <<  J(x + alpha_opt * dx)<< std::endl; 
            x.head(n_dofs * m) = f0;
            x.tail(n_dofs * m) = g0;
        
            // update
            if(J(x + dx) > J(x)){
              //std::cout << "entro: " << J(x + dx) << " " << J(x) << std::endl;
              double alpha_opt = backtracking(x, dx);
              f = f0 + alpha_opt * df;
              g = g0 + alpha_opt * dg;
              //std::cout << "alpha_opt: " << alpha_opt << std::endl;
            }else{
              f = f0 + df;
              g = g0 + dg;
            }

            f0 = f;
            g0 = g;
            
            //std::cout << "df (range) " << df.minCoeff() << " " << df.maxCoeff() << std::endl;   
            //std::cout << "f (range) " << f.minCoeff() << " " << f.maxCoeff() << std::endl;
            
        
            sse_new = (Psi*f - observations_).array().square().mean();
            std::cout << "sse_new   " << sse_new << std::endl;
            loss_new = sse_new + lambda*(g.transpose()*(Mt_*g))(0,0);
            std::cout << "J_old     "<< loss_old << std::endl;
            std::cout << "J_new     " << loss_new <<  "\n" << std::endl;
            // check

            stop = std::sqrt( (df.transpose()*(Mt_*df))(0,0) + (dg.transpose()*(Mt_*dg))(0,0) ) < tol || 
                  std::abs( (loss_new - loss_old)/loss_old ) < tol;
            
            std::cout << "||dx||_2 " << std::sqrt( (df.transpose()*(Mt_*df))(0,0) + (dg.transpose()*(Mt_*dg))(0,0) ) << std::endl;
            std::cout << "(J_new - J_old)/J_old " << std::abs( (loss_new - loss_old)/loss_old ) << std::endl;

            //stop = std::abs( loss_new - loss_old ) / loss_old  < tol; 
            //std::cout << "stop (Manzoni)    " << (std::sqrt( (df.transpose()*(Mt_*df))(0,0) + (dg.transpose()*(Mt_*dg))(0,0) ) < tol_)<< std::endl;
            //std::cout << "stop (Functional) " << (std::abs( loss_new - loss_old ) / loss_old < tol_) << std::endl; 
            //std::cout << "stop (sse)        " << (std::abs( sse_new - sse_old ) / sse_old < tol_ )<< std::endl;
            sse_old = sse_new;
            loss_old = loss_new;
            ++iter;
        }
        
        return std::make_pair(f, iter);
    };

    double J(const vector_t&x) const {
      vector_t f = x.head(n_dofs * m);
      vector_t g = x.tail(n_dofs * m);
      
      double SSE = (Psi*f - observations_).array().square().mean();
      double mis = lambda_*(g.transpose()*(Mt_*g))(0,0);
      return SSE + mis;
    }

    double backtracking(const vector_t& x, const vector_t& dx) const{

      double alpha_opt = alpha;
      vector_t f = x.head(n_dofs * m);
      vector_t g = x.tail(n_dofs * m);
      
      vector_t df = dx.head(n_dofs * m);
      vector_t dg = dx.tail(n_dofs * m);
      
      while( this->J(x + alpha_opt*dx) > (1.-c*alpha_opt)*this->J(x)){
        alpha_opt = beta * alpha_opt;
      }
      return alpha_opt;
    }

    void set_max_iter(int iters){ max_iter = iters;}
    void set_f_initial_guess(const vector_t& F0){f0 = F0;}
    void set_g_initial_guess(const vector_t& G0){g0 = G0;}
  // attribute
  const Triangulation<2, 2> &domain_;

  double mu_, r_;
  sparse_matrix_t A_; // [A]_{ij} = \int_D mu * grad(Psi_j)*grad(Psi_i) - alpha
                       // * Psi_j * Psi_i
  sparse_matrix_t M_; // Mass
  sparse_matrix_t At_;
  sparse_matrix_t psi;
  sparse_matrix_t Psi;
  sparse_matrix_t Mt_; // Spatio-Temporal Mass matrix
  sparse_matrix_t Psi_t_;

  sparse_matrix_t D, L;
  std::vector<sparse_matrix_t> Mw;
  std::vector<sparse_matrix_t> Gw;
  
  int n_dofs = 0;
  int n = 0;
  int m = 0;
  double deltaT = 0.;

  const matrix_t &locations_;
  const vector_t &time_locations_;
  vector_t observations_;

  double lambda_ = 1.;
  
  vector_t f0;
  vector_t g0;
        
  vector_t y_;
  vector_t u_;
  vector_t p_;

  vector_t IC;
  vector_t u_init_;
  double tol_ = 1e-4;
  int n_iter_;
  int max_iter = 10;

  double alpha = 1;
  double beta = 0.5;
  double c = 1e-4;

};



int main(int argc, char *argv[]){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "input/";
    std::string mesh_dir = "input/mesh/";
    std::string n_locs = std::string(argv[1]); // "100", "250" "500" "1000" 
    std::string sim = std::string(argv[2]); // 0, ..., 29
    std::string sim_dir = data_dir + n_locs  + "/" + sim + "/";

    std::string output_dir = "output-ic-esatta-newton/" +  n_locs  + "/" + sim + "/";
    std::string command_str = "mkdir -p " + output_dir;
    system(command_str.c_str());

    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    matrix_t test_locs = read_mtx<double>(mesh_dir + "test_locs.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);

    FeSpace Vh(unit_square, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);
    ZeroField<2> f;

    matrix_t time_mesh = read_mtx<double>(mesh_dir + "time_mesh.mtx");
    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    int n_times = time_mesh.rows();
    Triangulation<1, 1> T(time_mesh);

    matrix_t locs = read_mtx<double>(sim_dir + "locs.mtx");
    matrix_t obs = read_mtx<double>(sim_dir + "obs.mtx");
    
    matrix_t exact = read_mtx<double>(data_dir + "fisher_kpp.mtx"); 
    vector_t IC = exact.col(0);

    double mu = 0.01; // 0.001
    
    int n_dofs = nodes.rows();
    vector_t g_init = vector_t::Random(n_dofs*time_mesh.size()); // + vector_t::Ones(n_dofs*time_locs.size()));
    auto a = integral(unit_square)(mu * dot(grad(u), grad(v))); 
    auto F = integral(unit_square)(f * v);

    vector_t lambdas = vector_t::Ones(17);
    lambdas << 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 1e1, 5e1, 1e2, 5e2, 1e3;
    
    //vector_t lambdas = vector_t::Ones(4);
    //lambdas << 5, 10, 100, 1000; //1, 1e1;

    // vector_t lambdas = vector_t::Ones(1);
    // lambdas << 100.;

    vector_t reactions = vector_t::Ones(9);
    reactions << 0., 0.1, 0.25, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0;

    //vector_t reactions = vector_t::Ones(3);
    //reactions << 0.9, 1.0, 1.1;


    // ---------------------------------------------------------------------------------------------------------------------------
    Eigen::Matrix<double, Dynamic,2> grid = expand_grid( std::vector<vector_t>{lambdas, reactions} );

    std::cout << grid.transpose() << std::endl;
    
    Eigen::saveMarket(grid, output_dir + "cv_grids.mtx");

    int K = 10;
    KFoldCV cv_(locs, obs, K);
    vector_t cv_error = vector_t::Ones(K);

    auto mass = integral(unit_square)(u*v).assemble();
    sparse_matrix_t Im(n_times-1, n_times-1);
    Im.setIdentity();
    sparse_matrix_t Mt = kronecker(Im,mass);
    

  
    ScalarField<2> SSE;
    SSE = [&unit_square, &cv_, &time_mesh, &IC, &Vh, &mu, &n_times, &n_dofs, &g_init, &Im, &T, &a, &F, &time_locs](const vector_t& grid) -> double {
        double lambda = grid[0];
        double alpha = grid[1];
        if (lambda <= 0.0) return std::numeric_limits<double>::max();
        
        double res = 0;
        for(int k=0; k < cv_.k(); ++k){

            matrix_t locs = cv_.X_train(k);
            matrix_t obs = cv_.Y_train(k);
            vector_t  result = vector_t::Zero((n_times-1)*n_dofs);
            if( alpha > 0){

                vector_t response = obs.rightCols(n_times - 1).reshaped();
                GeoFrame data(unit_square, T);
                auto &l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{locs, time_locs});
                l.load_vec("y", response);
                SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));
                parabolic.fit(std::vector<double>{lambda, 1.0});

                // penso il newton vada implementato come un "ottimizzatore"
                // se esponesse un grad_old, update ecc potremmo usare backtracking linesearch?
                auto model = newton_method(unit_square, mu, alpha, locs, time_locs,
                                            response);
                model.set_f_initial_guess(parabolic.f());
                model.set_g_initial_guess(parabolic.misfit());
                model.set_state_initial_condition(IC);
                //model.set_control_initial_guess(g_init);
                model.set_max_iter(500);                
      
                auto res = model.solve(lambda); // BacktrackingLineSearch()

                matrix_t test_locs = cv_.X_test(k);
                matrix_t test_obs = cv_.Y_test(k).rightCols(n_times-1);
                sparse_matrix_t psi_test = internals::point_basis_eval(Vh, test_locs);
             
                sparse_matrix_t Psi_test = kronecker(Im,psi_test);
            
                result = std::get<0>(res); //model.y().tail((n_times-1)*n_dofs);
            
            }else{

                vector_t response = obs.rightCols(n_times - 1).reshaped();
                GeoFrame data(unit_square, T);
                auto &l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{locs, time_locs});
                l.load_vec("y", response);
                SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));
                parabolic.fit(std::vector<double>{lambda, 1.0});
                result = parabolic.f();
            }

            matrix_t test_locs = cv_.X_test(k);
            matrix_t test_obs = cv_.Y_test(k).rightCols(n_times-1);
            sparse_matrix_t psi_test = internals::point_basis_eval(Vh, test_locs);
            
            sparse_matrix_t Psi_test = kronecker(Im,psi_test);
            vector_t test_vals = Psi_test * result; //;matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
            res += std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());
        }
        
        std::cout << "lambda " << lambda << " alpha " << alpha <<  " cv_error = " << 1./cv_.k() * res<<  std::endl;
        return 1./cv_.k() * res;

    };

    GridSearch<2> optimizer;
    optimizer.optimize(SSE, grid); 
    
    //
    vector_t rmse = vector_t::Ones(1);
    matrix_t test_obs = read_mtx<double>(mesh_dir + "test_obs.mtx").rightCols(n_times-1);  
    auto psi_test = internals::point_basis_eval(Vh, test_locs);
    sparse_matrix_t Psi_test = kronecker(Im, psi_test);

    vector_t values = vector_t::Zero(optimizer.values().size());
    for(int i = 0; i < optimizer.values().size(); ++i) values[i] = optimizer.values()[i];
    
    Eigen::saveMarket(values, output_dir + "cv_errors.mtx");

    Eigen::Matrix<double, 2, 1> optimum(optimizer.optimum()[0], optimizer.optimum()[1]);

    Eigen::saveMarket(optimum, output_dir + "cv_optim.mtx");

    // PARA
    vector_t values_para = values.head(lambdas.size());
    Eigen::Index minRow;
    values_para.minCoeff(&minRow);
    double lambda_para = lambdas[minRow]; // !!!!! c'era messo "values" ma ci va "lambdas"
    std::cout << "opt (para) " << lambda_para << std::endl;
    std::cout << "values (para) " << values_para << std::endl;
    // ---- output kFold 

    vector_t response = obs.rightCols(n_times - 1).reshaped();
    GeoFrame data(unit_square, T);
    auto &l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{locs, time_locs});
    l.load_vec("y", response);
    SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));

    parabolic.fit( std::vector<double>{ lambda_para, 1.0} );
    vector_t test_vals = Psi_test * parabolic.f();
    rmse[0] = std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());

    std::cout << "rmse (para)" << rmse << std::endl;
    Eigen::saveMarket(rmse, output_dir + "rmse_diff_kfold.mtx");
    
    vector_t tmp = vector_t::Zero(n_dofs*n_times);
    tmp.head(n_dofs) = IC;
    tmp.tail(n_dofs*(n_times-1)) = parabolic.f();

    Eigen::saveMarket(tmp, output_dir + "estimate_diff_kfold.mtx");

    // NEWTOn

    std::cout << "opt " << optimizer.optimum() << std::endl;    
    std::cout << "values " << optimizer.values() << std::endl;
    
    if(optimum[1] > 0.){
      auto model = newton_method(unit_square, mu, optimum[1], locs, time_mesh,
                             response);
      model.set_f_initial_guess(parabolic.f());
      model.set_g_initial_guess(parabolic.misfit());
      model.set_state_initial_condition(IC);
      model.set_max_iter(500);
      auto res = model.solve(optimum[0]); // BacktrackingLineSearch()
      vector_t result = std::get<0>(res);
    
      test_vals = Psi_test * result; //;matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
      rmse[0] = std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());

      tmp.tail((n_times-1)*n_dofs) = result;
      std::cout << "rmse (NL)" << rmse << std::endl;
      Eigen::saveMarket(rmse, output_dir + "rmse_iterative.mtx");
      Eigen::saveMarket(tmp, output_dir + "estimate_iterative.mtx");
    }else{
      Eigen::saveMarket(rmse, output_dir + "rmse_iterative.mtx");
      Eigen::saveMarket(tmp, output_dir + "estimate_iterative.mtx");
    }
    return 0;
}