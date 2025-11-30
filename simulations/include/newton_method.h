#ifndef __FE_LS_FISHER_KPP_LUMPED_H__
#define __FE_LS_FISHER_KPP_LUMPED_H__

namespace fdapde {

template< int N>
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
    
    n_dofs = Vh.n_dofs();
    sparse_matrix_t psi(n,n_dofs);
    vector_t measure_vec = vector_t::Ones(n);

    if( incidence_matrix_.size() != 0){
        auto tmp = internals::areal_basis_eval(a.trial_space(), incidence_matrix_);
        psi = std::get<0>(tmp);
        measure_vec = std::get<1>(tmp);

    }else{
        psi = internals::point_basis_eval(a.trial_space(), locations_);
    }
    
    y_ = vector_t::Zero(n_dofs * m);
    u_ = y_;
    p_ = y_;

    sparse_matrix_t Im(m, m);
    Im.setIdentity();

    Mt_ = kronecker(Im, M_);
    Mt_.makeCompressed();
    
    Psi = kronecker(Im, psi);
    
    //matrix_t measure = measure_vec.asDiagonal();
    sparse_matrix_t measure(n,n);
    for(int i = 0; i < n; ++i) measure.insert(i,i) =  measure_vec[i];

    Measure = kronecker(Im, measure); 
    // std::cout << "regions " << n << std::endl;
    // std::cout << "Measure " <<  Measure.rows() << " " << Measure.cols() << std::endl;
    PsiT_D = Psi.transpose() * Measure;
    PsiT_D.pruned();

    //Measure = kronecker(Im, measure);
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

    f0 = vector_t::Zero(n_dofs*m); // initial guess -> iter 0 == "linearized"
    g0 = vector_t::Zero(n_dofs*m);
    // std::cout << "PsiT_D " <<  PsiT_D.rows() << " " << PsiT_D.cols() << std::endl; 
    // std::cout << "Psi " <<  Psi.rows() << " " << Psi.cols() << std::endl; 
    // std::cout << "Mt " << Mt_.rows() << " " << Mt_.cols() << std::endl;
    // std::cout << "finito " << std::endl;
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

  
  // POINTWISE
  newton_method(const Triangulation<N, N> &domain, double mu, double r,
                 const matrix_t &locations, const vector_t &time_locations,
                const vector_t &observations)
      : domain_(domain), mu_(mu), r_(r), locations_(locations),
        time_locations_(time_locations),
        deltaT(time_locations[1] - time_locations[0]), observations_(observations),
        n(locations_.rows()), m(time_locations.size()) {
      init();
        }

  // AREALE 
  newton_method(const Triangulation<N, N> &domain, double mu, double r,
                 const BinaryMatrix<Dynamic> &incidence_matrix, const vector_t &time_locations,
                const vector_t &observations)
      : domain_(domain), mu_(mu), r_(r), incidence_matrix_(incidence_matrix),
        time_locations_(time_locations),
        deltaT(time_locations[1] - time_locations[0]), observations_(observations),
        n(incidence_matrix.rows()), m(time_locations.size()) {
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
        
        //auto PsiTPsi = Psi.transpose()*Measure*Psi;
        
        update_weighted_mass(f0, Mw); 
        update_weighted_mass(g0, Gw); 
        
        vector_t u_ = vector_t::Zero(n_dofs * m); //no forcing
        std::cout << (M_*IC).rows() << " " << (M_*IC).cols() << std::endl;
        u_.segment(0, n_dofs) += (1.0 / deltaT) * (M_*IC);
        double loss_old = std::numeric_limits<double>::max();
        double loss_new = 0;
        loss_old = J(x);
        //std::cout << "J: " << loss_old << std::endl;
        
        double sse_old = std::numeric_limits<double>::max();
        double sse_new = 0; 
        
        //std::cout << observations_.rows() << " " << observations_.cols() << std::endl;
        while( !stop && iter < max_iter){
            std::cout << "\n\t iter " << iter << std::endl;
            update_weighted_mass(f0, Mw); 
            update_weighted_mass(g0, Gw); 

            sparse_matrix_t A_tilde = blockDiag(Mw);
            //std::cout << A_tilde.rows() << " " << A_tilde.cols() << std::endl;
            //range(A_tilde);

            sparse_matrix_t G_tilde = blockDiag(Gw);
            //range(G_tilde);
            
            //std::cout << G_tilde.rows() << " " << G_tilde.cols() << std::endl;
            //std::cout << Psi.rows() << " " << Psi.cols() << std::endl;

            auto S_00 = 1./(n * m)*PsiT_D *Psi + 2*r_*lambda*G_tilde;
            sparse_matrix_t S_10 = (D  + At_ -r_*Mt_ + 2*r_*A_tilde);

            sparse_matrix_t S_01 = lambda*S_10.transpose();

            SparseBlockMatrix<double, 2, 2> S ( S_00, S_01,
                                                S_10, -Mt_);
            
            invS_.compute(S);
            
            // right hand side
            b = vector_t::Zero(2 * n_dofs * m);

            b.head(n_dofs * m) = 1./(n * m) * PsiT_D * observations_ - 1./(n * m)* PsiT_D * (Psi*f0) - lambda*S_01*g0 ;
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
  const Triangulation<N, N> &domain_;

  double mu_, r_;
  sparse_matrix_t A_; // [A]_{ij} = \int_D mu * grad(Psi_j)*grad(Psi_i) - alpha
                       // * Psi_j * Psi_i
  sparse_matrix_t M_; // Mass
  sparse_matrix_t At_;
  sparse_matrix_t psi;
  sparse_matrix_t Psi;
  sparse_matrix_t Mt_; // Spatio-Temporal Mass matrix
  sparse_matrix_t Psi_t_;
  sparse_matrix_t  PsiT_D; 
  sparse_matrix_t Measure;
  sparse_matrix_t D, L;
  std::vector<sparse_matrix_t> Mw;
  std::vector<sparse_matrix_t> Gw;
  BinaryMatrix<Dynamic> incidence_matrix_;
  matrix_t locations_;

  int n_dofs = 0;
  int n = 0;
  int m = 0;
  double deltaT = 0.;

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

}

#endif