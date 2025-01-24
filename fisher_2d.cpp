#include "utils.h"
#include<fdaPDE/drivers.h>

using namespace fdapde;

struct fe_fisher_kpp {
	
	private:
    
    void init(){
    	
    	FeSpace Vh(domain_, P1<1>);
   		TrialFunction uh(Vh); 
  		TestFunction  vh(Vh);
  		
  		auto a = integral(domain_)(mu_*dot(grad(uh), grad(vh)) - alpha_*uh*vh);
   		auto mass = integral(domain_)(uh*vh);
   		
   		A_ = a.assemble();
   		A_.makeCompressed();
   		
   		M_ = mass.assemble();
   		M_.makeCompressed();
   		
   		Psi_ = internals::point_basis_eval(a.trial_space(), locations_);
   		Psi_.makeCompressed();
   		
   		PsiTPsi_ = Psi_.transpose() * Psi_;
  		PsiTPsi_.makeCompressed();
  		
  		n_dofs_ = Vh.n_dofs();
  		
  		y_ = VectorType::Zero(n_dofs_*n_time_locs_);
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
    using DenseSolverType  = Eigen::PartialPivLU<MatrixType>;
    
    fe_fisher_kpp(const Triangulation<2,2>& domain, double mu, double alpha,
    			  const MatrixType& locations, const VectorType& time_locations, const VectorType& observations): domain_(domain), mu_(mu), alpha_(alpha), 
    			  locations_(locations), time_locations_(time_locations),dT_(time_locations[1] - time_locations[0]),
    			  observations_(observations), n_locs_(locations.rows()), n_time_locs_(time_locations.size())
    			  {
    			  	init();
    			  }
    
    void solve(){
    	u_ = opt_.optimize(*this, u_init_);
    	y_ = state(u_);
    	p_ = adjoint(y_);
    }
    
    double J(const VectorType& y, const VectorType& u) const {
    	
    	double SSE = 0;
        for (int t = 0; t < n_time_locs_; ++t) {
            SSE += (obs(t) - Psi_*get_t(y,t)).squaredNorm();
        }
        return 1./(n_locs_ * n_time_locs_)*SSE + lambda_ * u.squaredNorm();
    
    }
    
    VectorType state(const VectorType& u) const{
    	
    	FeSpace Vh(domain_, P1<1>);
   		TrialFunction uh(Vh); 
  		TestFunction  vh(Vh);
  		FeFunction y_old(Vh);
  		
  		auto reac = integral(domain_)(alpha_ * y_old * uh * vh);
  		
    	VectorType y = VectorType::Zero(n_dofs_ * n_time_locs_);
    	y.block(0,0, n_dofs_, 1) = y0_;
    	
    	
    	for(int t = 0; t < n_time_locs_ - 1; t++){
    		y_old = get_t(y,t); 
    		auto R_ = reac.assemble();
    		R_.makeCompressed();
   		
    		SparseMatrixType S = 1./dT_*M_ + A_ + R_;
			S.prune(1e-10);
    		S.makeCompressed();
   			
    		SparseSolverType lin_solver(S);
    		lin_solver.factorize(S);
    		
    		VectorType b = 1./dT_*M_*y_old.coeff() + M_* get_t(u,t+1);
    		
    		y.block( (t+1)*n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b); 
    	}	
    	
    	return y;
    }
    
    VectorType adjoint(const VectorType& y) const{
    	
    	FeSpace Vh(domain_, P1<1>);
   		TrialFunction uh(Vh); 
  		TestFunction  vh(Vh);
  		FeFunction y_old(Vh);
  		
  		auto reac = integral(domain_)(alpha_ * y_old * uh * vh);
  		
    	VectorType p = VectorType::Zero(n_dofs_ * n_time_locs_); // p(T) == 0
    	
    	for(int t = n_time_locs_ - 1; t > 0; t--){
    		y_old = get_t(y,t-1);
    		auto R_ = reac.assemble();
    		R_.makeCompressed();
   		
    		SparseMatrixType S = 1./dT_*M_ + A_ + 2.*R_;
    		S.prune(1e-10);
    		S.makeCompressed();
   			
   			SparseSolverType lin_solver(S);
    		lin_solver.factorize(S);
    		
    		VectorType b = 1./dT_*M_*get_t(p,t) - 1./(n_time_locs_*n_locs_)*PsiTPsi_* get_t(y,t-1) + 
    											  1./(n_time_locs_*n_locs_)*Psi_.transpose()*obs(t-1);
    		
    		p.block( (t-1)*n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b); 
    	}
    
    	return p;
    }
    
    double operator()(const VectorType& u) const{

    	auto y = state(u);
    	return J(y,u);
    }
    
    // derive() deve restituire una lambda !
    std::function<VectorType(const VectorType& u)> derive(){
    	return [this](const VectorType& u){
    		auto y = state(u);
    		auto p = adjoint(y);
    		auto grad_u = alpha_ * Mt_ * u - Mt_ * p;
    		return grad_u;
    	};
    }
      
    VectorType y() const { 
    	return y_; 
    }
    
   VectorType u() const { 
    	return u_;
    }
    
   VectorType p() const { 
    	return p_; 
    }
    
    VectorType get_t(const VectorType& v, int t) const { return v.block(n_dofs_ * t, 0, n_dofs_, 1); }
    
    VectorType y(int t) const { return get_t(y(), t);}
    VectorType u(int t) const { return get_t(u(), t);}
    VectorType p(int t) const { return get_t(p(), t);}
    VectorType obs(int t)const { return observations_.block(n_locs_ * t, 0, n_locs_, 1);}
    
    // setters
    void set_state_initial_condition(const VectorType& y0){
    	y0_ = y0;
    }
    
    void set_control_initial_guess(const VectorType& u0){
    	u_init_ = u0;
    }
   
   	// optimization algorithm custom stopping criterion
    template <typename OptimizerType> bool opt_stopping_criterion(OptimizerType& opt) {
    	
    	bool stop;
    	
    	VectorType y_old = state(opt.x_old);
    	double loss_old = J(opt.x_old, y_old);
    	
    	VectorType y_new = state(opt.x_new);
    	double loss_new = J(opt.x_new, y_new);
		
		stop = std::abs((loss_new - loss_old) / loss_old) < tol_;
		return stop;
    }
   
    // attribute 
    const Triangulation<2,2>& domain_;
    
    double mu_, alpha_;
    SparseMatrixType A_; // [A]_{ij} = \int_D mu * grad(Psi_j)*grad(Psi_i) - alpha * Psi_j * Psi_i
    SparseMatrixType M_; // Mass
    
    SparseMatrixType Psi_; 
  	SparseMatrixType PsiTPsi_;
  	SparseMatrixType Mt_; // Spatio-Temporal Mass matrix	
  	
    int n_dofs_ = 0;
    int n_locs_ = 0;
    int n_time_locs_ = 0;
    double dT_ = 0.;
    
    const MatrixType& locations_;
    const VectorType& time_locations_;
    const VectorType& observations_;
    
    double lambda_ = 1.;
   
    VectorType y_;
    VectorType u_;
    VectorType p_;
    
	VectorType y0_;
	VectorType u_init_;
	double tol_ = 1e-6;
	BFGS<Eigen::Dynamic, WolfeLineSearch> opt_ = BFGS<Eigen::Dynamic, WolfeLineSearch>(30, tol_, 0.01);
};

int main(){
   std::string datadir = "data/";
   std::string meshdir = datadir + "mesh/";
   constexpr int local_dim = 2;
   using PointT = Eigen::Matrix<double, local_dim, 1>;
   
   Eigen::MatrixXd nodes = read_mtx<double>(meshdir +"nodes.mtx");
   Eigen::MatrixXi elements = read_mtx<int>(meshdir +"elements.mtx");
   Eigen::MatrixXi boundary = read_mtx<int>(meshdir +"boundary.mtx");
   Triangulation<local_dim, local_dim> unit_square = Triangulation<2, 2>(nodes, elements, boundary);
   
   
   std::cout << nodes.rows() << " " << nodes.cols() << std::endl;
   
   double mu = 0.1;
   double alpha = 3.0;
   
   int n_observations = read_txt<int>(datadir + "n_locations.txt")(1,0);
   int nsim = 30;
   int sim = 0;
   
   Eigen::MatrixXd locations = read_txt<double>(datadir + std::to_string(n_observations) + "/" + std::to_string(sim)  + "/" + "locs.txt"); 
   Eigen::MatrixXd obsMat = read_txt<double>(datadir + std::to_string(n_observations) + "/" + std::to_string(sim) + "/" + "obs.txt");
   std::cout << locations.rows() << " " << locations.cols() << std::endl;
   
   Eigen::VectorXd observations = Eigen::Map<Eigen::VectorXd>(obsMat.data(), obsMat.size());
   std::cout << observations.rows() << " " << observations.cols() << std::endl;
   
   Eigen::VectorXd time_locations = read_txt<double>(datadir + "time_locations.txt");
   
   std::cout << time_locations.rows() << " " << time_locations.cols() << std::endl;
   
   auto model = fe_fisher_kpp(unit_square, mu, alpha, locations, time_locations, observations);
   
   
   Eigen::VectorXd IC = read_txt<double>(datadir + "exact.txt").col(0);
   std::cout << IC.rows() << " " << IC.cols() << std::endl;
   
   model.set_state_initial_condition(IC);
   
   Eigen::VectorXd u0 = read_txt<double>(datadir + "u_guess_rand.txt");
   model.set_control_initial_guess(u0);
   
   model.solve();
   
   std::cout << model.y().rows() <<  " " << model.y().cols() << std::endl; 
   std::cout <<"n iterations: " <<  model.opt_.n_iter() << std::endl;
   
   std::cout <<"state: " << model.y().minCoeff() <<  " " << model.y().maxCoeff() << std::endl; 
   std::cout <<"control: " << model.u().minCoeff() <<  " " << model.u().maxCoeff() << std::endl; 
   std::cout <<"adjoint: " << model.p().minCoeff() <<  " " << model.p().maxCoeff() << std::endl; 
   Eigen::saveMarket(model.y(), datadir + std::to_string(n_observations) +"/" + "y.mtx");   
   return 0;
}

