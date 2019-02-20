
class PairStyle {
  public:
    int type1, type2 ;
    double *u, **f, **vir ;
    double *rho1, *rho2 ;
    complex<double> *u_k, **f_k, **vir_k ;

    PairStyle( int, double*, double* ) ; 
    ~PairStyle() ;
    
    double calc_energy() ;
} ;
