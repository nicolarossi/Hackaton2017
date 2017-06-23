// BS_hackaton.cpp : Defines the entry point for the console application.
#include <omp.h>
#include <iostream>
#include <array>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <chrono>
#include <vector>

#define INITIAL_DV_STEP (1e-4)
#define MAX_THREAD (232)

#define M_PI 3.14159265358979323846
//#define M_SQRT1_2 0.707106781186547524401


/* To optimize the call_price function, note: it's 4 are mutually exclusive */
//#define FUNCTION
//#define USING_MACRO
//#define REORDER_1
#define REORDER_2
//#define USING_TEMPLATE

/* A strategy to parallelize and speedup the single stock loop in a simple way define or not to enable/disable it. 
  disabling this and enabling DYNAMIC_DV you test loop also with a dv dynamic (disabling both test with dv static
*/
#define INDIRECT_INDEX
#define DYNAMIC_DV


/* A strategie to "increase" the efficience of parallel loop */
//#define USING_REDUCTION



/* DON'T TOUCH UNDER THIS ROW */

#if defined FUNCTION
double norm_cdf(double value) {
  return 0.5 * erfc(-value * M_SQRT1_2);
}

double d1(double S, double K, double r, double v, double T)
{
  return (log(S / K) + r*T) / (v*(pow(T, 0.5))) + 0.5*v*(pow(T, 0.5));
}

double d2(double S, double K, double r, double v, double T)
{
  return d1(S, K, r, v, T) - v*(pow(T, 0.5));
}
#elif defined  REORDER_1
 #define norm_cdf( value) ( 0.5 * erfc(-value * M_SQRT1_2))
 #define d1( S,  K,  r,  v,  T,v_pow) (( (log(S / K) + r*T) / (v_pow)) + 0.5*v_pow )
 #define d2( S,  K,  r,  v, T,v_pow)  (( (log(S / K) + r*T) / (v_pow)) - 0.5*v_pow )
#elif defined REORDER_2
 #define norm_cdf( value) ( 0.5 * erfc(-value * M_SQRT1_2))
 #define d1( S,  K,  r,  v,  T,v_pow,log_S_K_p_rT) (( (log_S_K_p_rT) / (v_pow)) + 0.5*v_pow )
 #define d2( S,  K,  r,  v, T,v_pow,log_S_K_p_rT)  (( (log_S_K_p_rT) / (v_pow)) - 0.5*v_pow )
#elif defined USING_MACRO
  #define norm_cdf( value) ( 0.5 * erfc(-value * M_SQRT1_2))
  #define d1( S,  K,  r,  v,  T) ( (log(S / K) + r*T) / (v*(pow(T, 0.5))) + 0.5*v*(pow(T, 0.5)))
  #define d2( S,  K,  r,  v, T)  ( d1(S, K, r, v, T) - v*(pow(T, 0.5)))
#elif defined USING_TEMPLATE
template<typename T>
inline T norm_cdf( T value) {
  return ( 0.5 * erfc(-value * M_SQRT1_2));
}

template<typename T>
inline T d1( T S,  T K,  T r,  T v,  T t) {
  return ( (log(S / K) + r*t) / (v*(pow(t, 0.5))) + 0.5*v*(pow(t, 0.5)));
}

template<typename T>
inline T d2( T S, T  K, T r,  T v, T t)  {
  return ( d1(S, K, r, v, t) - v*(pow(t, 0.5)));
}
#endif

int n_max;

double exp_rT, v_pow, log_S_K_p_RT,c_pow;

// vanilla european call price based on
// underlying S, strike K, risk-free rate r, volatility of
// underlying sigma and time to maturity T
double _call_price(double S, double K, double r, double v, double T) 
{
#if defined REORDER_1
  double v_pow=v*(pow(T, 0.5));
  return S * norm_cdf(d1(S, K, r, v, T,v_pow)) - K*exp(-r*T) * norm_cdf(d2(S, K, r, v, T,v_pow));
#elif defined REORDER_2
  extern double c_pow,log_S_K_p_RT,exp_rT;
  double v_pow=v*c_pow;
  //.....
  //  pow_T=(pow(T, 0.5) 
  //  double v_pow=v*pow_T
  //  double log_S_K_p_RT=(log(S/K)+r*T)*0.5;
  //  double exp_rT=exp(-r*T);
  
  return S * norm_cdf(d1(S, K, r, v, T,v_pow,log_S_K_p_RT)) - K*exp_rT * norm_cdf(d2(S, K, r, v, T,v_pow,log_S_K_p_RT));
#else
  return S * norm_cdf(d1(S, K, r, v, T)) - K*exp(-r*T) * norm_cdf(d2(S, K, r, v, T));
#endif
}

// vanilla european put price based on
// underlying S, strike K, risk-free rate r, volatility of
// underlying sigma and time to maturity T
double put_price(double S, double K, double r, double v, double T) 
{
#ifdef REORDER_1
  double v_pow=v*(pow(T, 0.5));
  return -S * norm_cdf(-d1(S, K, r, v, T,v_pow)) + K*exp(-r*T) * norm_cdf(-d2(S, K, r, v, T,v_pow));
#elif defined REORDER_2
  std::cerr<< "NOT YET IMPLEMENTED at line "<<__LINE__<<"\n";
  exit(-1);
#else
  return -S * norm_cdf(-d1(S, K, r, v, T)) + K*exp(-r*T) * norm_cdf(-d2(S, K, r, v, T));
#endif
}

//---- We look for call_to_find from a,b with step dv

class solution_t {
public :
  double v,call;
  static double value_to_find;

  //---
  //solution_t(double d): v(-1.),call(std::numeric_limits<double>::max()/10.),value_to_find(d) {
  /*  solution_t(double to_find): v(-1.),call(-1.),value_to_find(to_find) {
  }*/

  solution_t(): v(-1.),call(std::numeric_limits<double>::max()) {
  }

  //---
  void call_price(double S, double K, double r, double v, double T);
  
  solution_t& operator+( solution_t &rhs);
  bool operator<(solution_t &rhs) ;

  friend std::ostream& operator<<(std::ostream&os,solution_t &obj){
    os<<" v: \t"<<obj.v<<" \t call: \t"<<obj.call<<" \t value_to_find: \t"<<obj.value_to_find;
    return os;
  }
  
};  

double  solution_t::value_to_find;

void solution_t::call_price(double S, double K, double r, double v, double T) {
  this->call=_call_price( S,  K,  r,  v,  T);
  this->v=v;
}

bool solution_t::operator<(solution_t &rhs) {
  //--- so we return the max volatility that match
  //  std::cerr << "\n\n Asking if \n" <<*this<<" \nis less then \n"<<rhs<<"\n";
  
  if (fabs(this->call-this->value_to_find)==0.) return true; //--- a < b  iff ||a-x||= 0. || or 3||a-x||<||b-x||
  return (fabs(this->call-solution_t::value_to_find) < fabs(rhs.call-solution_t::value_to_find));
}


//---

struct step_t {
  double v;
  double dv;
  
  step_t(double _v,double _dv):v(_v),dv(_dv){};  
};

std::vector<step_t> step;

void  __create_partitioned_loop(double a,double b,double dv) {
  double v;
  //--- Calculare number of step
  for (v=a; v<=b ; v+=dv) {
    step.emplace_back(v,dv);
    
    //---
    dv=( v*INITIAL_DV_STEP) ;
    dv=(dv<INITIAL_DV_STEP?INITIAL_DV_STEP:dv); /*---- In the initial phase of the loop the cycle is increased of a FIXED_STEP (e.g 1e-3 
						  after that phase the step is relative to the v eg when v = 1000 dv= 1 (v*dv)
						*/
  }  
}

solution_t find_volatility(double call_to_find,double S,double K, double r, double T, double a, double b,double dv) {//, solution_t best_finded) {
  solution_t tmp_rv;
  solution_t::value_to_find=call_to_find;
  
  double v=a;
  
#ifdef INDIRECT_INDEX

#ifdef  USING_REDUCTION
    solution_t rv_min;
    
#pragma omp declare reduction ( best_solution : solution_t : omp_out = omp_in < omp_out ? omp_in : omp_out ) 
#pragma omp parallel for private(tmp_rv,v,dv) reduction(best_solution : rv_min)
  for (int i=0;i<step.size();i++) {
    v=(step[i].v);
    tmp_rv.call_price(S, K,  r,  v,  T);

    if (tmp_rv<rv_min) {
      rv_min=tmp_rv;
    }
    dv=step[i].dv;
  }
  
  return rv_min;

#else
  std::vector<solution_t> rv(MAX_THREAD);
  
  // for (v=a+dv; v<=b ; v+=dv) {
  //  for (auto it=step.begin();it!=step.end();++it) {
#pragma omp parallel for shared(rv),private(tmp_rv,v,dv)
  for (int i=0;i<step.size();i++) {
    v=(step[i].v);
    tmp_rv.call_price(S, K,  r,  v,  T);
    int who_i_am=omp_get_thread_num();
    //    std::cerr<<" debug "<<i<<"\t"<< tmp_rv.call<<"\t\t";
    if (tmp_rv<rv[who_i_am]) {
      rv[who_i_am]=tmp_rv;
      //      std::cerr<<" cambiamo ";
    }
    //    std::cerr<<" \n";
    dv=step[i].dv;
  }
  
  int index_min=0;
  for (int i=0;i<omp_get_max_threads();i++){
    if (rv[i]<rv[index_min]){
      index_min=i;
    }
  }
  
  return rv[index_min];
#endif

#else
#ifdef DYNAMIC_DV
  solution_t rv;

  for (v=a; v<=b ; v+=dv){
    tmp_rv.call_price(S, K,  r,  v,  T);
    if (tmp_rv<rv) {
      rv=tmp_rv;
    }
    
    dv=( v*INITIAL_DV_STEP) ;
    dv=(dv<INITIAL_DV_STEP?INITIAL_DV_STEP:dv); /*---- In the initial phase of the loop the cycle is increased of a FIXED_STEP (e.g 1e-3 
						  after that phase the step is relative to the v   eg when v = 1000 dv= 1 (v*dv)
						*/
  }
  
  return rv;
#else
  solution_t rv;

  for (v=a; v<=b ; v+=dv){
    tmp_rv.call_price(S, K,  r,  v,  T);
    if (tmp_rv<rv) {
      rv=tmp_rv;
    }
    
    dv=( INITIAL_DV_STEP) ;
  }
  
  return rv;
#endif  
#endif
  
}

int main(int argc, char **argv) 
{
  std::cerr<<std::setprecision(7);
  std::vector<double> time_result;
  int readed=0;
  int processed=0;
  
  if (argc<2) {
    std::cerr<<" Usage\n\t"<<argv[0]<<" <number of threads> [n_sample]\n";
    exit(-1);
  }
  int nt=atoi(argv[1]);
  omp_set_num_threads(nt);
  n_max=omp_get_max_threads();
  std::cerr << " n_max = "<<n_max<<"\n";
  
  std::ifstream fi;
  std::istringstream is;
  std::string line;
  std::array<std::string, 20> columns;
  std::vector<std::array<std::string,20>> sample;
  int n_sample=(argc>2?atoi(argv[2]):100);
  
  double sum_call=0.;
  

  fi.open("Sample_SPX_20151001_to_20151030_WITH_SPXW.csv");
  while (getline(fi, line, '\n')) {
    if ( n_sample <=readed ) {
      break;
    }
    readed++;
    
    is.str(line);
    for (size_t index = 0; getline(is, columns[index], ','); index++);
    is.clear();
    
    sample.push_back(columns);
  }

  
  auto startTime = std::chrono::high_resolution_clock::now();
  
  double a,b,dv;
  
  a=0.; b=150.;
  dv=INITIAL_DV_STEP;

#ifdef INDIRECT_INDEX
  __create_partitioned_loop(a,b,dv);
#endif

  double r = 0.02; // risk-free interest rate (2%)
  double v_real = 0.2; // volatility of the underlying (20%)
  double T = 10./252.; // 1 year to expiry
  exp_rT=exp(r*T);
  c_pow=pow(T,0.5);
  log_S_K_p_RT;
  for (auto it=sample.begin();it!=sample.end();++it) {
    if (processed%10==0) std::cerr<< " Processing row ="<< processed <<"\n";

    double S = atof( (*it)[1 /* underlying_last */].c_str());
    double K = atof((*it)[8 /* strike */].c_str());
    
    processed++;
    
    double call_real = _call_price(S, K, r, v_real, T);

    //--- 
    
    solution_t best_finded;//(call_real);
    auto internal_startTime = std::chrono::high_resolution_clock::now();

    log_S_K_p_RT=log(S/K)+r*T;
    best_finded=find_volatility(call_real, S, K,  r,  T,  a,b, dv);
    
    auto internal_endTime = std::chrono::high_resolution_clock::now();
    
    
    //    if (0){
    if (processed%10==1) {
      std::cerr<< " \nCall to find : "<< call_real << "";
      std::cerr<< " Best value of v found : "<< best_finded.v << "";
      std::cerr<< " Call founded : "<< best_finded.call << "\n";  
    }
    //      sum_call+=best_finded.call;
    
    if (0) {
      std::cerr<< " ERROR WE DIDN'T FIND THE VOLATILITY FOR row ="<< readed <<"\n";
      exit(-1);
    } else {
      double t = std::chrono::duration_cast<std::chrono::milliseconds>(internal_endTime-internal_startTime).count();
      
      
      time_result.push_back(t);
      //      std::cerr<< " Time : "<< t << " msec\n";
    }
    //    std::cerr<< "\n";
    
  }

  
  fi.close();


   // Exercise 2: try to find 'v' (volatility) by reversing the call_price function.
  //  Say you start with the value for 'call', S, K, r, T, but not v. And try to 
    //  solve what the input value 'v' should be to get the value for 'call'. 
    //
    // This means you need to invoke call_price() many times (trying different
    //  values for volatility) until you find a return value that matches the
    //  call value.

  auto endTime = std::chrono::high_resolution_clock::now();
  std::cerr << " Readed "<< readed << " processed "<< processed<<"\n";
  auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count();
 
   std::cerr << " Time :"<<diff<<"\n";
   std::cerr << " Horse Power :"<<((double)processed/(double)diff)<<"\n";

   auto sum=time_result[0];
   sum=0;
   for (auto it=time_result.begin();it!=time_result.end();++it){
     sum+=*it;
   }
   auto avg=(sum/time_result.size());
   std::cerr << " Avg "<< avg<<" msec\n";

   sum=0;
   for (auto it=time_result.begin();it!=time_result.end();++it){
     sum+=(*it-avg)*(*it-avg);
   }

   std::cerr << " Sigma "<< sqrt(sum/time_result.size())<<"\n";

   

    return 0;
}
