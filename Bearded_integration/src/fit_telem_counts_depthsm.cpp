
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function to make sparse SPDE precision matrix
// Inputs :
// logkappa : log ( kappa ) parameter value
// logtau : log ( tau ) parameter value
// M0 , M1 , M2: these sparse matrices are output from :
// R:: INLA :: inla . spde2 . matern () $ param . inla $M*
template <class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0 ,
                              SparseMatrix <Type> M1 , SparseMatrix <Type> M2) {
  SparseMatrix <Type> Q;
  Type kappa2 = exp (2. * logkappa);
  Type kappa4 = kappa2 * kappa2;
  Q = pow (exp(logtau), 2.) * (kappa4*M0 + Type(2.0)*kappa2*M1+M2);
  return Q;
}

// helper function to use the same penalized complexity prior on
// matern params that is used in INLA

template <class Type>
Type dPCPriSPDE (Type logtau, Type logkappa,
                  Type matern_par_a, Type matern_par_b,
                  Type matern_par_c, Type matern_par_d,
                  int give_log=0)
{

  Type penalty; // prior contribution to jnll
  
  Type d=2.; // dimension
  Type lambda1= -log(matern_par_b)*pow(matern_par_a, d/2.);
  Type lambda2= -log(matern_par_d)/matern_par_c;
  Type range = sqrt(8.0)/exp(logkappa);
  Type sigma = 1.0/sqrt(4.0*3.14159265359*exp(2.0*logtau)*
    exp(2.0*logkappa));
  
  penalty=(-d/2. - 1.)*log(range)-lambda1*pow(range,-d/2.)-lambda2*sigma;
  // Note : (rho , sigma ) --> (x= log kappa , y= log tau ) -->
  // transforms : rho = sqrt (8) /e^x & sigma = 1/( sqrt (4 pi)*e^x*e^y)
  // --> Jacobian : |J| propto e^( -y -2x)
  Type jacobian = -logtau-2.0*logkappa;
  penalty += jacobian;
  
  if( give_log ) return penalty; else return exp(penalty);
}



//dmvnorm function with pre-calculated matrix inverse, log determinant - ignore normalizing constant
template <class Type>
Type neg_log_dmvnorm(vector<Type> x, vector<Type> mu, matrix<Type> omega) {
 vector<Type> ss = x-mu;
 return ((0.5*ss*(omega*ss)).sum());
}

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}


template<class Type>
Type objective_function<Type>::operator() ()
{


  // Data
  DATA_MATRIX( Y_i );       	// (n_s * n_ind) MATRIX of telemetry location counts for each data type at each sampled location 
                            //NB: this will have to redone in more realistic situations where each surface has a different # of samples
  DATA_MATRIX( X );   // fixed effects design matrix 
  DATA_MATRIX( X_pred );   // design matrix for predictions
  DATA_VECTOR( Offset);
  DATA_IVECTOR(Season); //season for each observation
  
  DATA_SPARSE_MATRIX( S );
  DATA_MATRIX(X_sm);
  DATA_MATRIX(X_sm_pred);
  DATA_MATRIX(depthReport);
  
  // normalization flag - used for speed -up
  DATA_INTEGER(flag); // flag == 0 => no data contribution added to jnll
  
  // Indices
  DATA_INTEGER( n_ind ); // Number of unique individual/season combos
  DATA_INTEGER( n_s ); // Number of grid cells
  DATA_INTEGER( n_yrs); //number of years for predictions
  DATA_INTEGER( n_mesh ); //mesh points in space mesh
  DATA_IVECTOR( N_samp); //number of observations of each individual


  // SPDE objects  - assume these are all the same for each surface
  DATA_SPARSE_MATRIX ( M0 );
  DATA_SPARSE_MATRIX ( M1 );
  DATA_SPARSE_MATRIX ( M2 );
  DATA_SPARSE_MATRIX ( A );
  //DATA_IVECTOR ( Mesh_index ); //which mesh points are associated with grid cell centroids
  
  // Options
  DATA_VECTOR ( options );
  // options [0] == 1 : use normalization trick

  // Prior specifications
  DATA_VECTOR( matern_pri );  //NOT CURRENTLY USED
  // matern_pri = c(a, b, c, d): P( range < a) = b; P( sigma > c) = d
  Type matern_par_a = matern_pri[0]; // range limit : rho0
  Type matern_par_b = matern_pri[1]; // range prob : alpha_rho
  Type matern_par_c = matern_pri[2]; // field sd limit : sigma0
  Type matern_par_d = matern_pri[3]; // field sd prob : alpha_sigma

 
  // Parameters 
  PARAMETER_VECTOR( Beta );              // fixed effects on UD counts
  PARAMETER_VECTOR( Beta_sm);         //smooth depth effects
  PARAMETER( log_tau_eta );      //spatial precision for overall UD
  PARAMETER( log_kappa_eta );
  PARAMETER( log_tau_xi );      //spatial precision for individual tracks
  PARAMETER( log_kappa_xi );
  PARAMETER( log_lambda ); //penalty on smooth effects

  PARAMETER_MATRIX( Eta_s );           // seasonal spatial mean REs
  PARAMETER_MATRIX( Xi_s );           // spatial individual REs 

    // derived sizes
  //int n_b = X_s.row(0).size();
  
  // global stuff
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  Type jnll = 0.0;
  int n_seasons = Eta_s.row(0).size();
  int n_b_sm = S.row(0).size();
  Type lambda=exp(log_lambda);

  // Make spatial precision matrix
  SparseMatrix <Type > Q_eta = spde_Q( log_kappa_eta , log_tau_eta , M0 , M1 , M2);
  SparseMatrix <Type > Q_xi;
  // the random effects . we do this first so to do the
  // normalization outside of every optimization step
  if( options [0] == 1){
    // then we are not calculating the normalizing constant in the inner opt
    // that norm constant means taking an expensive determinant of Q_ss
    Q_xi = spde_Q( log_kappa_xi , log_tau_xi , M0 , M1 , M2);
    for(int iind=0;iind<n_ind;iind++){
      jnll_comp(0) += GMRF(Q_xi , false )(Xi_s.col(iind));
    }
    for(int iseas=0;iseas<n_seasons;iseas++){
      jnll_comp(0) += GMRF(Q_eta , false )(Eta_s.col(iseas));
    }
    // return without data ll contrib to avoid unneccesary log ( det (Q)) calcs
    if ( flag == 0) return jnll_comp(0) ; //see documentation for TMB::normalize
  }
  else {
    Q_xi = spde_Q( log_kappa_xi, log_tau_xi , M0 , M1 , M2);
    for(int iind=0;iind<n_ind;iind++){
      jnll_comp(0) += GMRF(Q_xi)(Xi_s.col(iind));
    }
    for(int iseas=0;iseas<n_seasons;iseas++){
      jnll_comp(0) += GMRF(Q_eta)(Eta_s.col(iseas));
    }
  }
  
  matrix<Type> Eta_pred(n_s,n_seasons);
  matrix<Type> Xi_pred(n_s,n_ind);
  for(int iind=0;iind<n_ind;iind++){
    Xi_pred.col(iind)=A*Xi_s.col(iind);
  }
  for(int iseas=0;iseas<n_seasons;iseas++){
    Eta_pred.col(iseas)=A*Eta_s.col(iseas);
  }

  int n_tot = n_ind * n_s;
  vector<Type> Lin_pred = X * Beta + X_sm * Beta_sm;
  matrix<Type> Lambda(n_s,n_ind);
  int cur_pl = 0;
  for(int iind=0; iind<n_ind;iind++){
    for(int is=0;is<n_s;is++){
      Lambda(is,iind) = exp(Offset(iind)+Lin_pred(cur_pl+is)+ Xi_pred(is,iind) + Eta_pred(is,Season(iind)));
      jnll_comp(1)+= Lambda(is,iind)-Y_i(is,iind)*log(Lambda(is,iind));   //pois model
    }
    cur_pl+=n_s;
  }

  vector<Type> A_Eta(n_s);
  vector<Type> XB = X_pred * Beta + X_sm_pred *Beta_sm;
  vector<Type> XB_tmp(n_s);
  
  matrix<Type> log_UD_s(n_s,n_yrs*n_seasons);
  for(int iyr=0; iyr<n_yrs; iyr++){
    for(int iseas=0;iseas<n_seasons;iseas++){
      A_Eta= A*Eta_s.col(iseas);
      XB_tmp = XB.segment(iyr*(n_s*n_seasons)+iseas*n_s,n_s);
      log_UD_s.col(iyr*n_seasons+iseas)=XB_tmp + A_Eta;
    }
  }
  
  // //priors
  // //jnll_comp(2) -= dnorm (Beta , beta_pri [0] , beta_pri [1] , true );

  jnll_comp(2) -= Type(0.5)*n_b_sm*log_lambda - 0.5*lambda*GMRF(S).Quadform(Beta_sm); //note from Devin: m_i would need to be rank(S) if S_i not full rank (e.g. thin plate)
  jnll_comp(2) -= Type(-0.95)*log_lambda-0.005*lambda;  //gamma(0.05,0.005) prior like in Jagam

  vector<Type> splineForReport = depthReport*Beta_sm;
  
  // // Total objective
  jnll = jnll_comp.sum();

  //std::cout<<"jnll_comp "<<jnll_comp<<"\n";
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( Xi_s );
  REPORT( Eta_s );
  REPORT( Xi_pred);
  REPORT( Eta_pred);
  REPORT( log_tau_eta );
  REPORT( log_kappa_eta );
  REPORT( log_tau_xi );
  REPORT( log_kappa_xi );
  REPORT( lambda);
  REPORT( Lambda );
  REPORT( Beta );
  REPORT( Beta_sm);
  REPORT( log_UD_s);

  ADREPORT(log_UD_s);
  ADREPORT(splineForReport);

  return jnll;
}
