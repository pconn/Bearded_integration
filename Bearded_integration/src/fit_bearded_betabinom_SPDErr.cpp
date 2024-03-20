#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

//dmvnorm function with pre-calculated matrix inverse omitting normalizing constant
template <class Type>
Type neg_log_dmvnorm(vector<Type> x, vector<Type> mu, matrix<Type> omega) {
  vector<Type> mu_diff = x-mu;
  return ((0.5*mu_diff*(omega*mu_diff)).sum());
}

// binomial log likelihood omitting normalizing constant and trapping numerical underflows
template<class Type>
Type dbinom_kern_log(Type n, Type x, Type p){
  Type p_num = p;
  if(p==0) p_num = 0.00000001;
  if(p==1) p_num = 0.99999999;
  //return x*log(p+p0)+(n-x)*log(1.0-p+p1);
  return x*log(p_num)+(n-x)*log(1.0-p_num);
}

// binomial log likelihood omitting normalizing constant and trapping numerical underflows
// this version uses n, x/n instead of n, x
template<class Type>
Type dbinom_kern_log2(Type n, Type p_obs, Type p){
  Type p_num = p;
  if(p==0) p_num = 0.00000001;
  if(p==1) p_num = 0.99999999;
  return n*p_obs*log(p_num)+n*(1-p_obs)*log(1.0-p_num);
}

// log beta function
template<class Type>
Type lbeta(Type z1, Type z2){
  return lgamma(z1)+lgamma(z2)-lgamma(z1+z2);
}

// log beta function
template<class Type>
Type dbetabinom_kern_log(Type n, Type y, Type mu, Type rho){
  if(rho<0.00000001) rho = 0.00000001;
  if(rho>0.99999999) rho = 0.99999999;
  if(mu<0.00000001) mu = 0.00000001;
  if(mu>0.99999999) mu = 0.99999999;
  Type tmp = (1.0-rho)/rho;
  Type alpha = mu*tmp;
  Type beta = (1-mu)*tmp;
  return lbeta(y+alpha,n-y+beta)-lbeta(alpha,beta);
}

// poisson approx to binom
template<class Type>
Type dpois_kern_log(Type n, Type x, Type p){
  Type lambda = n*p;
  return x*log(lambda)-lambda;
}

// poisson log likelihood kernel
template<class Type>
Type dpois_kern_log2(Type x, Type lambda){
  return x*log(lambda)-lambda;
}

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

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  // Data + fixed inputs
  DATA_MATRIX( Acc_n );      // holds number of acoustic duty cycle intervals (binomial trials) per "season" (1 record per duty cycle type)
  DATA_MATRIX( Acc_k );    // holds number of acoustic duty cycle intervals for which calls were recorded (binomial successes)
  DATA_MATRIX( Area_acc );  // each row gives the proportion of each grid cell for which seals are "detectable" by the mooring represented by each detection record
  DATA_VECTOR( Duty_lengths);  //length of each duty cycle in seconds
  DATA_IVECTOR( Acc_tstep);  // holds time step for each acoustic detection summary
  DATA_MATRIX( X_acc );           // design matrix to explain variation in acoustic detections        
  DATA_SCALAR( log_ckmr_N );   // CKMR estimate (log scale)
  DATA_SCALAR( log_ckmr_se );  // standard error of CKMR estimate (log scale)
  DATA_MATRIX(UD_mean_adj);  // utilization distribution, standardized by divided by "mean_UD"
  DATA_SCALAR(mean_UD);  // for matching Pi_st with standardized UD
  DATA_IVECTOR(which_tsteps_UD); //don't fit UDs from years with no data
  DATA_MATRIX(W_st); // inverse variance for UD observations - used in weighted sum-of-squares
  DATA_IVECTOR( which_Bering );  //which cells of larger grid are part of the Bering Sea
  DATA_IVECTOR( which_Chukchi );  //which cells of larger grid are part of the Chukchi Sea
  DATA_IVECTOR( which_Beaufort );  //which cells of larger grid are part of the Chukchi Sea
  DATA_VECTOR( log_2012_N );  //log of abundance values predicted for each Bering grid cell in 2012
  DATA_VECTOR( log_2013_N );  //log of abundance values predicted for each Bering grid cell in 2013
  DATA_VECTOR( log_2016_N );  //log of abundance values predicted for each Chukchi grid cell in 2016
  DATA_VECTOR( log_2021_N );  //log of abundance values predicted for each Chukchi grid cell in 2016
  DATA_VECTOR( log_2012_se );
  DATA_VECTOR( log_2013_se );
  DATA_VECTOR( log_2016_se );
  DATA_VECTOR( log_2021_se );
  DATA_VECTOR( N_2012 );  //abundance values predicted for each Bering grid cell in 2012 (for chi-sq)
  DATA_VECTOR( N_2013 );  //abundance values predicted for each Bering grid cell in 2012 (for chi-sq)
  DATA_VECTOR( N_2016 );  //abundance values predicted for each Bering grid cell in 2012 (for chi-sq)
  DATA_VECTOR( N_2021 );  //abundance values predicted for each Bering grid cell in 2012 (for chi-sq)
  DATA_ARRAY( X_ice );  //each [,,i] matrix is a design matrix for time step i
  DATA_VECTOR( Land_cover ); //proportion of each cell that is composed of land
  DATA_SPARSE_MATRIX( designMatrixForReport );//Design matrix for report of splines
  DATA_SPARSE_MATRIX( S );//Penalization block diagonal matrix 
  DATA_IVECTOR( Sdims );   //Dimensions of each smooth
  DATA_MATRIX( X_s );  //design matrix associated with spline effects
  DATA_VECTOR( Wts );  //likelihood weights for composite likelihood. Order: UD, CKMR, Aerial, Acoustics
  DATA_INTEGER( n_tsteps);  //number of time steps (# number of years * # of seasons)
  DATA_INTEGER( n_seasons);  //number of seasons per year
  DATA_INTEGER(est_acc); //only add in acoustic likelihood if est_acc=1
  // SPDE objects
  DATA_INTEGER( n_eta ); //mesh points in space mesh
  DATA_SPARSE_MATRIX ( M0 );
  DATA_SPARSE_MATRIX ( M1 );
  DATA_SPARSE_MATRIX ( M2 );
  DATA_SPARSE_MATRIX ( A ); //for predictive process model
  // normalization flag - used for speed -up
  DATA_INTEGER(flag); // flag == 0 => no data contribution added to jnll
  
  

  // Parameters 
  PARAMETER( log_N );
  PARAMETER_VECTOR ( Beta_s );  //parameters governing spline effects of abundance intensity
  PARAMETER_VECTOR ( Beta_ice );  //parameters governing quadratic sea ice effects
  PARAMETER_VECTOR ( Beta_acc );  //parameters governing acoustic detection process
  PARAMETER_VECTOR( log_lambda );      //penalization parameters for spline model 
  PARAMETER( log_tau );      //spatial precision
  PARAMETER( log_kappa );
  PARAMETER_VECTOR( Etainput_s );           // spatial random effects
  PARAMETER( logit_rho );  //overdispersion parameter for mu / rho parameterization of beta binom

  // derived sizes
  int n_b = X_s.row(0).size();
  //int n_b_ice = X_ice.col(0).row(0).size();
  int n_lambda = log_lambda.size();
  int n_cells = X_s.col(0).size();
  //vector<Type> ice_dim = X_ice.dim;
  //int n_b_ice = ice_dim(1);
  vector<int> dim_ice = X_ice.dim;
  int n_b_ice = X_ice.dim(1);
  int n_chukchi = which_Chukchi.size();
  int n_bering = which_Bering.size();
  int n_beaufort = which_Beaufort.size();
  int n_beta_acc = X_acc.row(0).size();
  int n_acc_det = Acc_n.col(0).size(); 
  int n_duty = Duty_lengths.size();
  int n_tsteps_UD = which_tsteps_UD.size();
  int n_yrs = n_tsteps/n_seasons;

  // global stuff
  vector<Type> jnll_comp(6);  //neg. log likelihood components: UD, CKMR, Aerial, Acoustics, spline prior
  jnll_comp.setZero();
  Type jnll = 0.0;
  Type N = exp(log_N);
  vector<Type> lambda(n_lambda);
  Type rho = 1.0/(1.0+exp(-logit_rho));
  for(int ilam=0;ilam<n_lambda;ilam++)lambda(ilam)= exp(log_lambda(ilam));
  
  // Make spatial precision matrix
  SparseMatrix <Type > Q_ss = spde_Q( log_kappa , log_tau , M0 , M1 , M2);

  // Transform some of our parameters
  Type sp_range = sqrt(8.0) / exp(log_kappa);
  Type sp_sigma = 1.0 / sqrt (4.0 * 3.14159265359 * exp(2.0*log_tau) * exp(2.0*log_kappa));

  // the random effects . we do this first so to do the
  // normalization outside of every optimization step
    // not calculating the normalizing constant in the inner opt
    // that norm constant means taking an expensive determinant of Q_ss
  int start_pl=0;
  for(int iseas=0;iseas<n_seasons;iseas++){
    jnll_comp(5) += GMRF(Q_ss , false )(Etainput_s.segment(start_pl,n_eta));
    start_pl = start_pl+n_eta;
  }

    // return without data ll contrib to avoid unneccesary log ( det (Q)) calcs
  if ( flag == 0) return jnll_comp(5);

  //calculate expected abundance
  matrix<Type> Pi_st(n_cells,n_tsteps);  //proportion in each cell - this will be matched to UDs
  matrix<Type> Z_st(n_cells,n_tsteps);
  matrix<Type> logZ_st(n_cells,n_tsteps);
  matrix<Type> logPi_st(n_cells,n_tsteps);
  vector<Type> fixed_eff_s = X_s * Beta_s;  // time-invariant spline linear predictor
  vector<Type> linpredZ_st(n_cells);

  Type cur_sum;
  Type small = 0.00000001;
  int it=0;
  vector<Type> RE(n_cells);
  vector<Type> Eta_s(n_eta);

  for(int iseas=0;iseas<n_seasons;iseas++){
    start_pl = iseas*n_eta;
    Eta_s = Etainput_s.segment(start_pl,n_eta);
    RE = A * Eta_s;
    for(int iyr = 0; iyr<n_yrs; iyr++){
      it = iyr*2+iseas;
      for(int is=0; is<n_cells; is++){
        linpredZ_st(is) = fixed_eff_s(is)+RE(is);
        for(int ib = 0; ib<n_b_ice; ib++){
          linpredZ_st(is) += X_ice(is,ib,it)*Beta_ice(ib);
        }
      }
      cur_sum = 0;
      for(int is=0; is<n_cells;is++){
        Pi_st(is,it) = (1-Land_cover(is))*exp( linpredZ_st(is) ) + small;
        cur_sum+=Pi_st(is,it);
      }
      for(int is=0; is<n_cells;is++){
        Pi_st(is,it) = Pi_st(is,it)/cur_sum;
        Z_st(is,it) = N*Pi_st(is,it);
        logPi_st(is,it) = log(Pi_st(is,it));
        logZ_st(is,it) = log_N + logPi_st(is,it);
      }
    }
  }

  //UD observation process
  matrix<Type> Pi_st_mean_adj = Pi_st/mean_UD;
  for(int it=0; it<n_tsteps_UD; it++){
    for(int is=0; is<n_cells; is++){
      jnll_comp(0)+=W_st(is,it)*pow(UD_mean_adj(is,which_tsteps_UD(it))-Pi_st_mean_adj(is,which_tsteps_UD(it)),2.0);  // inv-var weighted sum of squares...
    }
  }
  jnll_comp(0)=Wts(0)*jnll_comp(0);

  //ckmr observation process
  //jnll_comp(1) = -Wts(1)*dnorm(log_N,log_ckmr_N,log_ckmr_se,1);  //log likelihood - normal on log-scale
  jnll_comp(1) = -dnorm(log_N,log_ckmr_N,Wts(1)*log_ckmr_se,1);  //log likelihood - normal on log-scale
  //jnll_comp(1) = 0.0;

  //aerial survey counts
  vector<Type> log_E_N_2012(n_bering);
  vector<Type> log_E_N_2013(n_bering);
  vector<Type> log_E_N_2016(n_chukchi);
  vector<Type> log_E_N_2021(n_beaufort);
  vector<Type> Z_score_2012(n_bering);
  vector<Type> Z_score_2013(n_bering);
  vector<Type> Z_score_2016(n_chukchi);
  vector<Type> Z_score_2021(n_beaufort);
  Type chisq_aerial=0.0;

  using namespace density;

  int t_step1 = 16;  //corresponds to model time step associated with april/may 2012
  int t_step2 = 18;  //corresponds to model time step associated with april/may 2013
  int t_step3 = 24;  //corresponds to model time step associated with april/may 2016
  int t_step4 = 34;  //corresponds to model time step associated with april/may 2021
  //log_2012_se = log_2012_se*Wts(2);
  //log_2013_se = log_2013_se*Wts(2);
  //log_2016_se = log_2016_se*Wts(2);
  //log_2021_se = log_2021_se*Wts(2);

  for(int i=0; i<n_bering; i++){
    log_E_N_2012(i)=logZ_st(which_Bering(i),t_step1);
    log_E_N_2013(i)=logZ_st(which_Bering(i),t_step2);
    chisq_aerial += pow(N_2012(i)-Z_st(which_Bering(i),t_step1),2.0)/Z_st(which_Bering(i),t_step1);
    chisq_aerial += pow(N_2013(i)-Z_st(which_Bering(i),t_step2),2.0)/Z_st(which_Bering(i),t_step2);
    Z_score_2012(i) = (log_E_N_2012(i)-log_2012_N(i))/log_2012_se(i);
    Z_score_2013(i) = (log_E_N_2013(i)-log_2013_N(i))/log_2013_se(i);
    jnll_comp(2) -= dnorm(log_E_N_2012(i),log_2012_N(i),log_2012_se(i),1);  //speed up w/ kernel
    jnll_comp(2) -= dnorm(log_E_N_2013(i),log_2013_N(i),log_2013_se(i),1);
  }
  for(int i=0; i<n_chukchi; i++){
    log_E_N_2016(i)=logZ_st(which_Chukchi(i),t_step3);
    Z_score_2016(i) = (log_E_N_2016(i)-log_2016_N(i))/log_2016_se(i);
    jnll_comp(2) -= dnorm(log_E_N_2016(i),log_2016_N(i),log_2016_se(i),1);
    chisq_aerial += pow(N_2016(i)-Z_st(which_Chukchi(i),t_step3),2.0)/Z_st(which_Chukchi(i),t_step3);
  }
  for(int i=0; i<n_beaufort; i++){
   log_E_N_2021(i)=logZ_st(which_Beaufort(i),t_step4);
   Z_score_2021(i) = (log_E_N_2021(i)-log_2021_N(i))/log_2021_se(i);
   jnll_comp(2) -= dnorm(log_E_N_2021(i),log_2021_N(i),log_2021_se(i),1);
   chisq_aerial += pow(N_2021(i)-Z_st(which_Beaufort(i),t_step3),2.0)/Z_st(which_Beaufort(i),t_step3);
  }
  
  jnll_comp(2)=jnll_comp(2)*Wts(2);

  //acoustics
  Type Expect;
  vector<Type> Rate = exp(X_acc * Beta_acc);
  vector<Type> CurZ(n_cells);
  vector<Type> Moor_prop(n_cells);
  vector<Type> N_acc(n_acc_det);
  matrix<Type> Prob_acc(n_acc_det,n_duty);
  matrix<Type> Acc_p(n_acc_det,n_duty);
  matrix<Type> Acc_neff(n_acc_det,n_duty);
  Acc_p.setZero();
  Type chisq_acc = 0.0;
  matrix<Type> Chisq_acc(n_acc_det,n_duty);
  Chisq_acc.setZero();
  if(est_acc==1){
    for(int i=0; i<n_acc_det; i++){
      CurZ = Z_st.col(Acc_tstep(i));
      Moor_prop = Area_acc.row(i);
      N_acc(i) = (CurZ*Moor_prop).sum();  //number of "detectable" seals
      for(int iduty = 0; iduty<n_duty; iduty++){  //5 different duty cycle types
        if(Acc_n(i,iduty)>0){
          Acc_p(i,iduty) = Acc_k(i,iduty)/Acc_n(i,iduty);
          Prob_acc(i,iduty) = 1-exp(-Rate(i)*N_acc(i)*Duty_lengths(iduty));
          //Acc_neff(i,iduty) = Acc_n(i,iduty)*Wts(3);
          Chisq_acc(i,iduty) = Acc_n(i,iduty)*pow(Acc_p(i,iduty)-Prob_acc(i,iduty),2.0)/Prob_acc(i,iduty);
          chisq_acc += Acc_n(i,iduty)*pow(Acc_p(i,iduty)-Prob_acc(i,iduty),2.0)/Prob_acc(i,iduty);
          //jnll_comp(3) += dpois_kern_log(Acc_n(i,iduty),Acc_k(i,iduty),Prob_acc(i,iduty));
          //jnll_comp(3) += dbinom_kern_log2(Acc_neff(i,iduty),Acc_p(i,iduty),Prob_acc(i,iduty));
          jnll_comp(3) -= dbetabinom_kern_log(Acc_n(i,iduty),Acc_k(i,iduty),Prob_acc(i,iduty),rho);
        }
      }
    }
  }
  jnll_comp(3) = Wts(3)*jnll_comp(3);
  

  // //spline prior
  // int k=0;  // Counter
  // for(int i=0;i<n_lambda;i++){
  //   int m_i = Sdims(i);
  //   vector<Type> beta_i = Beta_s.segment(k,m_i);       // Recover betai
  //   SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
  //   jnll_comp(4) -= Type(0.5)*m_i*log_lambda(i) - 0.5*lambda(i)*GMRF(S_i).Quadform(beta_i); //note from Devin: m_i would need to be rank(S) if S_i not full rank (e.g. thin plate)
  //   k += m_i;
  //   jnll_comp(4) -= (Type(-0.95)*log_lambda(i)-0.005*exp(log_lambda(i)));  //gamma(0.05,0.005) prior like in Jagam
  // }
  // jnll_comp(4)=Wts(4)*jnll_comp(4);


  jnll = jnll_comp.sum();

  REPORT(jnll_comp);
  REPORT(Pi_st_mean_adj);
  REPORT(logPi_st);
  REPORT(Z_st);
  REPORT(N);
  REPORT(Beta_s);
  REPORT(Beta_ice);
  REPORT(lambda);
  REPORT(fixed_eff_s);
  REPORT(linpredZ_st);
  REPORT(Pi_st);
  REPORT(logZ_st);
  REPORT(log_N);
  REPORT(N);
  REPORT(log_E_N_2012);
  REPORT(log_E_N_2013);
  REPORT(log_E_N_2016);
  REPORT(log_E_N_2021);
  REPORT(N_acc);
  REPORT(Prob_acc);
  REPORT(Beta_acc);
  REPORT(chisq_aerial);
  REPORT(chisq_acc);
  REPORT(Chisq_acc);
  REPORT(Z_score_2012);
  REPORT(Z_score_2013);
  REPORT(Z_score_2016);
  REPORT(Z_score_2021);
  REPORT(RE);
  REPORT( log_tau );
  REPORT( log_kappa );
  REPORT( sp_sigma );
  REPORT( sp_range );
  REPORT( Etainput_s );
  REPORT( rho );
 
  //std::cout<<"Range "<<Range_eta<<"\n";
  //ADREPORT(log_Z_s);
  //ADREPORT(N);
  //ADREPORT(Beta_s);
  //ADREPORT(Beta_ice);
  //ADREPORT(Beta_acc);
  
  return jnll;
}
