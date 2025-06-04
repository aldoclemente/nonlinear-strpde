#include "../utils.h"
#include <fdaPDE/fdapde.h>

using namespace fdapde;

/*
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
    // I_ *= dT_; // ???

    Mt_ = kronecker(I_, M_);
    Mt_.makeCompressed();

    // Node
    std::cout << " init \n " << std::endl;
    std::cout << time_locations_ << std::endl;
    std::cout << observations_.rows() << std::endl;
    std::cout << n_time_locs_ << std::endl;
    std::cout << n_locs_ << std::endl;
  }

public:
  using VectorType = Eigen::Matrix<double, Dynamic, 1>;
  using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
  using SparseMatrixType = Eigen::SparseMatrix<double>;
  using DiagonalMatrixType = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
  using SparseSolverType = Eigen::SparseLU<SparseMatrixType>;
  using DenseSolverType = Eigen::PartialPivLU<MatrixType>;
  using PointT = Eigen::Matrix<double, 2, 1>;

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
    // beta_ = lambda_ * n_locs_ * n_time_locs_;
    u_ = opt_.optimize(*this, u_init_);
    y_ = state(u_);
    p_ = adjoint(y_);
  }

  double J(const VectorType &y, const VectorType &u) const {
    double SSE = 0;
    for (int t = 0; t < n_time_locs_; ++t) {
      SSE += (obs(t) - Psi_ * get_t(y, t)).squaredNorm();
    }
    return 1. / (n_locs_ * n_time_locs_) * SSE + lambda_ * u.squaredNorm();
    // return SSE + beta_ * u.squaredNorm();
  }

  VectorType state(const VectorType &u) const {
    // std::cout << "\t solving state " << std::endl;
    FeSpace Vh(domain_, P1<1>);
    auto &dof_handler = Vh.dof_handler();
    dof_handler.set_dirichlet_constraint(4, g_1);

    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    VectorType y = VectorType::Zero(n_dofs_ * n_time_locs_);
    y.block(0, 0, n_dofs_, 1) = y0_;

    for (int t = 0; t < n_time_locs_ - 1; t++) {
      y_old = get_t(y, t);
      auto R_ = reac.assemble();
      R_.makeCompressed();

      SparseMatrixType S = 1. / dT_ * M_ + A_ + R_; //+ eps_;
      dof_handler.enforce_constraints(S);
      S.makeCompressed();

      SparseSolverType lin_solver(S);
      lin_solver.factorize(S);

      VectorType b = 1. / dT_ * M_ * y_old.coeff() + M_ * get_t(u, t + 1);
      dof_handler.enforce_constraints(b);

      // std::cout << "rhs:\n" << b << std::endl;
      y.block((t + 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return y;
  }

  VectorType adjoint(const VectorType &y) const {
    // std::cout << "\t solving adjoint " << std::endl;
    FeSpace Vh(domain_, P1<1>);
    auto &dof_handler = Vh.dof_handler();
    dof_handler.set_dirichlet_constraint(4 g_0);

    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    VectorType p = VectorType::Zero(n_dofs_ * n_time_locs_); // p(T) == 0

    for (int t = n_time_locs_ - 1; t > 0; t--) {
      y_old = get_t(y, t - 1);
      // y_old = get_t(y, t);
      auto R_ = reac.assemble();
      R_.makeCompressed();

      SparseMatrixType S = 1. / dT_ * M_ + A_ + 2. * R_; // + eps_;
      dof_handler.enforce_constraints(S);
      S.makeCompressed();

      SparseSolverType lin_solver(S);
      lin_solver.factorize(S);

      VectorType b =
          1. / dT_ * M_ * get_t(p, t) -
          1. / (n_time_locs_ * n_locs_) * PsiTPsi_ * get_t(y, t - 1) +
          1. / (n_time_locs_ * n_locs_) * Psi_.transpose() * obs(t - 1);
      // VectorType b = 1. / dT_ * M_ * get_t(p, t) - PsiTPsi_ * get_t(y, t -
      // 1)+
      //                  Psi_.transpose() * obs(t - 1);
      //  dof_handler.enforce_constraints(S,b);
      dof_handler.enforce_constraints(b);
      p.block((t - 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return p;
  }

  double operator()(const VectorType &u) const {

    auto y = state(u);
    return J(y, u);
  }

  // derive() deve restituire una lambda !
  std::function<VectorType(const VectorType &u)> derive() {
    return [this](const VectorType &u) {
      auto y = state(u);
      auto p = adjoint(y);
      auto grad_u = lambda_ * Mt_ * u - Mt_ * p;

      // auto grad_u = u - p;
      std::cout << "\t update gradient " << std::endl;
      // auto grad_u = beta_ * Mt_ * u - Mt_ * p;
      return grad_u;
    };
  }

  VectorType y() const { return y_; }

  VectorType u() const { return u_; }

  VectorType p() const { return p_; }

  VectorType get_t(const VectorType &v, int t) const {
    return v.block(n_dofs_ * t, 0, n_dofs_, 1);
  }

  VectorType y(int t) const { return get_t(y(), t); }
  VectorType u(int t) const { return get_t(u(), t); }
  VectorType p(int t) const { return get_t(p(), t); }
  VectorType obs(int t) const {
    return observations_.block(n_locs_ * t, 0, n_locs_, 1);
  }

  // setters
  void set_state_initial_condition(const VectorType &y0) { y0_ = y0; }

  void set_control_initial_guess(const VectorType &u0) { u_init_ = u0; }

  // optimization algorithm custom stopping criterion

  // STOP IF !!!
  template <typename OptimizerType> bool stop_if(OptimizerType &opt) {

    bool stop;

    VectorType y_old = state(opt.x_old);
    // double loss_old = J(opt.x_old, y_old);

    VectorType y_new = state(opt.x_new);
    // double loss_new = J(opt.x_new, y_new);

    // if ((y_old - y_new).squaredNorm() / y_old.squaredNorm() > tol_)
    //   return false;

    double SSE_old = 0;
    double SSE_new = 0;
    for (int t = 0; t < n_time_locs_; ++t) {
      SSE_new += (obs(t) - Psi_ * get_t(y_new, t)).squaredNorm();
      SSE_old += (obs(t) - Psi_ * get_t(y_old, t)).squaredNorm();
    }

    double PEN_old = opt.x_old.squaredNorm();
    double PEN_new = opt.x_new.squaredNorm();

    double loss_old =
        1. / (n_locs_ * n_time_locs_) * SSE_old + lambda_ * PEN_old;
    double loss_new =
        1. / (n_locs_ * n_time_locs_) * SSE_new + lambda_ * PEN_new;

    double dloss = std::abs((loss_new - loss_old) / loss_old);
    stop = dloss < tol_;
    return stop;


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

  ScalarField<2, decltype([](const PointT &p) { return 0; })> g_0;
  ScalarField<2, decltype([](const PointT &p) { return 1.; })> g_1;

  double mu_, alpha_;
  SparseMatrixType A_; // [A]_{ij} = \int_D mu * grad(Psi_j)*grad(Psi_i) - alpha
                       // * Psi_j * Psi_i
  SparseMatrixType M_; // Mass

  SparseMatrixType Psi_;
  SparseMatrixType PsiTPsi_;
  SparseMatrixType Mt_; // Spatio-Temporal Mass matrix

  int n_dofs_ = 0;

  int n_locs_ = 0;
  int n_time_locs_ = 0;
  double dT_ = 0.;

  const MatrixType &locations_;
  const VectorType &time_locations_;
  const VectorType &observations_;

  double lambda_ = 1.;
  // double beta_ = lambda_ * n_locs_ * n_time_locs_;
  VectorType y_;
  VectorType u_;
  VectorType p_;

  VectorType y0_;
  VectorType u_init_;
  double tol_ = 1e-5;
  // con BFGS non va bene inizializzare u_init_ = 0 obv
  // BFGS<Eigen::Dynamic, BacktrackingLineSearch> opt_ =
  //    BFGS<Eigen::Dynamic, BacktrackingLineSearch>(500, tol_, 0.01);
  GradientDescent<Eigen::Dynamic, BacktrackingLineSearch> opt_ =
      GradientDescent<Eigen::Dynamic, BacktrackingLineSearch>(1000, tol_, 0.01);
  //   GradientDescent<Eigen::Dynamic> opt_ =
  //      GradientDescent<Eigen::Dynamic>(5000, tol_, 1e-2);
};
*/

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
    u_ = opt_.optimize(*this, u_init_);
    y_ = state(u_);
    p_ = adjoint(y_);
  }

  double J(const VectorType &y, const VectorType &u) const {

    double SSE = 0;
    for (int t = 0; t < n_time_locs_; ++t) {
      SSE += (obs(t) - Psi_ * get_t(y, t)).squaredNorm();
    }
    return 1. / (n_locs_ * n_time_locs_) * SSE + lambda_ * u.squaredNorm();
  }

  VectorType state(const VectorType &u) const {

    FeSpace Vh(domain_, P1<1>);
    auto &dof_handler = Vh.dof_handler();
    // Vh.impose_dirichlet_constraint(4, g_0);
    //  dof_handler.set_dirichlet_constraint(BoundaryAll, g_0);
    dof_handler.set_dirichlet_constraint(4, g_1);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    VectorType y = VectorType::Zero(n_dofs_ * n_time_locs_);
    y.block(0, 0, n_dofs_, 1) = y0_;

    for (int t = 0; t < n_time_locs_ - 1; t++) {
      y_old = get_t(y, t);
      auto R_ = reac.assemble();
      R_.makeCompressed();

      SparseMatrixType S = 1. / dT_ * M_ + A_ + R_;
      S.prune(1e-10);
      dof_handler.enforce_constraints(S);
      S.makeCompressed();

      SparseSolverType lin_solver(S);
      lin_solver.factorize(S);

      VectorType b = 1. / dT_ * M_ * y_old.coeff() + M_ * get_t(u, t + 1);
      dof_handler.enforce_constraints(b);
      y.block((t + 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return y;
  }

  VectorType adjoint(const VectorType &y) const {

    FeSpace Vh(domain_, P1<1>);
    auto &dof_handler = Vh.dof_handler();
    // Vh.impose_dirichlet_constraint(4, g_0); // questo solo sui modelli statis
    //  dof_handler.set_dirichlet_constraint(BoundaryAll, g_0);
    dof_handler.set_dirichlet_constraint(4, g_0);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    VectorType p = VectorType::Zero(n_dofs_ * n_time_locs_); // p(T) == 0

    for (int t = n_time_locs_ - 1; t > 0; t--) {
      y_old = get_t(y, t - 1);
      auto R_ = reac.assemble();
      R_.makeCompressed();

      SparseMatrixType S = 1. / dT_ * M_ + A_ + 2. * R_;
      S.prune(1e-10);
      dof_handler.enforce_constraints(S);
      S.makeCompressed();

      SparseSolverType lin_solver(S);
      lin_solver.factorize(S);

      VectorType b =
          1. / dT_ * M_ * get_t(p, t) -
          1. / (n_time_locs_ * n_locs_) * PsiTPsi_ * get_t(y, t - 1) +
          1. / (n_time_locs_ * n_locs_) * Psi_.transpose() * obs(t - 1);
      dof_handler.enforce_constraints(b);
      p.block((t - 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return p;
  }

  double operator()(const VectorType &u) const {

    auto y = state(u);
    return J(y, u);
  }

  // derive() deve restituire una lambda !
  std::function<VectorType(const VectorType &u)> derive() {
    return [this](const VectorType &u) {
      auto y = state(u);
      auto p = adjoint(y);
      auto grad_u = lambda_ * Mt_ * u - Mt_ * p; // qua ci va lambda e non alpha
      return grad_u;
    };
  }

  std::function<VectorType(const VectorType &u)> gradient() {
    return [this](const VectorType &u) {
      auto y = state(u);
      auto p = adjoint(y);
      auto grad_u = lambda_ * Mt_ * u - Mt_ * p; // qua ci va lambda e non alpha
      return grad_u;
    };
  }

  VectorType y() const { return y_; }

  VectorType u() const { return u_; }

  VectorType p() const { return p_; }

  VectorType get_t(const VectorType &v, int t) const {
    return v.block(n_dofs_ * t, 0, n_dofs_, 1);
  }

  VectorType y(int t) const { return get_t(y(), t); }
  VectorType u(int t) const { return get_t(u(), t); }
  VectorType p(int t) const { return get_t(p(), t); }
  VectorType obs(int t) const {
    return observations_.block(n_locs_ * t, 0, n_locs_, 1);
  }

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

  int n_dofs_ = 0;
  int n_locs_ = 0;
  int n_time_locs_ = 0;
  double dT_ = 0.;

  const MatrixType &locations_;
  const VectorType &time_locations_;
  const VectorType &observations_;

  double lambda_ = 1.;
  ScalarField<2,
              decltype([](const Eigen::Matrix<double, 2, 1> &p) { return 0; })>
      g_0;

  ScalarField<2,
              decltype([](const Eigen::Matrix<double, 2, 1> &p) { return 1; })>
      g_1;
  VectorType y_;
  VectorType u_;
  VectorType p_;

  VectorType y0_;
  VectorType u_init_;
  double tol_ = 1e-5;
  BFGS<Eigen::Dynamic, BacktrackingLineSearch> opt_ =
      BFGS<Eigen::Dynamic, BacktrackingLineSearch>(100, tol_, 0.01);
  // GradientDescent<Eigen::Dynamic, BacktrackingLineSearch> opt_ =
  //     GradientDescent<Eigen::Dynamic, BacktrackingLineSearch>(1000, tol_,
  //     0.01);
};

int main(int argc, char *argv[]) {
  std::cout << "running simulation 1" << std::endl;
  std::string datadir = "data/";
  std::string meshdir = datadir + "mesh/";
  constexpr int local_dim = 2;
  using PointT = Eigen::Matrix<double, local_dim, 1>;
  Eigen::MatrixXd nodes = read_mtx<double>(meshdir + "nodes.mtx");
  Eigen::MatrixXi elements = read_mtx<int>(meshdir + "elements.mtx");
  Eigen::MatrixXi boundary = read_mtx<int>(meshdir + "boundary.mtx");

  Triangulation<local_dim, local_dim> unit_square =
      Triangulation<2, 2>(nodes, elements, boundary);

  unit_square.mark_boundary(
      /* as = */ 4,
      /* where = */ [](const typename Triangulation<2, 2>::EdgeType &edge) {
        return (edge.node(0)[0] == 0.0 &&
                edge.node(1)[0] == 0.0); // left side ( x = 0 )
      });

  // ???

  // for(int i = 1; i < 5; ++i){
  /* for (auto it = unit_square.boundary_begin(4);
       it != unit_square.boundary_end(4); ++it) {
    std::cout << "(" << (*it).node(0)[0] << ", " << (*it).node(0)[1] << "); ("
              << (*it).node(1)[0] << ", " << (*it).node(1)[1] << ") \t"
              << (*it).marker() << std::endl;
    // all and only the boundary edges marked as
  }
   */
  //}

  double mu = 0.1;
  double alpha = 1.;

  int n_observations = read_TXT<int>(datadir + "n_locations.txt")(
      std::stoi(argv[1]), 0); // 0, ..., 4
  int nsim = 30;
  int sim = std::stoi(argv[2]); // 0, ..., 29

  Eigen::MatrixXd locations =
      read_TXT<double>(datadir + std::to_string(n_observations) + "/" +
                       std::to_string(sim) + "/" + "locs.txt");

  Eigen::MatrixXd obsMat =
      read_TXT<double>(datadir + std::to_string(n_observations) + "/" +
                       std::to_string(sim) + "/" + "obs.txt");
  std::cout << locations.rows() << " " << locations.cols() << std::endl;

  // Eigen::MatrixXd obsMat2 = obsMat.rightCols(obsMat.cols() - 1);
  Eigen::VectorXd observations =
      Eigen::Map<Eigen::VectorXd>(obsMat.data(), obsMat.size());
  std::cout << observations.rows() << " " << observations.cols() << std::endl;

  Eigen::VectorXd time_locations =
      read_TXT<double>(datadir + "time_locations.txt");

  std::cout << time_locations.rows() << " " << time_locations.cols()
            << std::endl;

  auto model = fe_fisher_kpp(unit_square, mu, alpha, locations, time_locations,
                             observations);

  // eh ...
  Eigen::VectorXd IC = read_TXT<double>(datadir + "exact.txt").col(0);
  std::cout << IC.rows() << " " << IC.cols() << std::endl;
  double lambda = 1;
  if (argc > 3) {
    lambda = std::stod(argv[3]);
  }

  // initialization with SRPDE -----------
  /*
  Triangulation<1, 1> T(time_locations.tail(time_locations.size() - 1));
  GeoFrame data(unit_square, T);
  auto &l = data.insert_scalar_layer<POINT, POINT>(
      "layer", std::pair{locations, MESH_NODES});
  l.load_vec("y",
             observations.tail((time_locations.size() - 1) * locations.rows()));
  std::cout << l << std::endl;
  FeSpace Vh(unit_square, P1<1>);
  TrialFunction u(Vh);
  TestFunction v(Vh);
  ZeroField<2> f;
  auto &dof_handler = Vh.dof_handler();
  auto a = integral(unit_square)(mu * dot(grad(u), grad(v)) - alpha * u * v);
  auto F = integral(unit_square)(f * v);
  SRPDE srpde("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));

  srpde.fit(lambda, lambda);
  Eigen::saveMarket(srpde.f(), datadir + std::to_string(n_observations) + "/" +
                                   std::to_string(sim) + "/" +
                                   "y_parabolic.mtx");
  Eigen::saveMarket(srpde.misfit(), datadir + std::to_string(n_observations) +
                                        "/" + std::to_string(sim) + "/" +
                                        "mis_parabolic.mtx");
  */
  // --- fine ---
  ScalarField<2, decltype([](const PointT &p) { return 1.; })> g_1;

  std::cout << IC.minCoeff() << " " << IC.maxCoeff() << std::endl;
  // dof_handler.set_dirichlet_constraint(/* on = */ 4, /* data = */ g_1);
  // dof_handler.enforce_constraints(IC);
  model.set_state_initial_condition(IC);

  Eigen::VectorXd u0 = read_TXT<double>(datadir + "u_guess_rand.txt");
  // model.set_control_initial_guess(srpde.misfit());
  model.set_control_initial_guess(u0);
  std::cout << "lambda = " << lambda << std::endl;
  // std::cout << "dim u0: " << u0.rows() << " " << u0.cols() << std::endl;
  //  Eigen::VectorXd U0 =
  //      Eigen::MatrixXd::Ones(nodes.rows() * time_locations.size(), 1);
  //  model.set_control_initial_guess(U0);
  model.solve(lambda);

  std::cout << model.y().rows() << " " << model.y().cols() << std::endl;
  std::cout << "n iterations: " << model.opt_.n_iter() << std::endl;

  std::cout << "state: " << model.y().minCoeff() << " " << model.y().maxCoeff()
            << std::endl;
  std::cout << "control: " << model.u().minCoeff() << " "
            << model.u().maxCoeff() << std::endl;
  std::cout << "adjoint: " << model.p().minCoeff() << " "
            << model.p().maxCoeff() << std::endl;

  vector2txt(std::vector<int>{model.opt_.n_iter()},
             datadir + std::to_string(n_observations) + "/" +
                 std::to_string(sim) + "/" + "n_iter.txt");
  Eigen::saveMarket(model.y(), datadir + std::to_string(n_observations) + "/" +
                                   std::to_string(sim) + "/" + "y.mtx");

  std::cout << "simulation 1 ended" << std::endl;
  return 0;
}
