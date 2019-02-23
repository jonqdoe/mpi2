#ifndef PAIR_STYLE
#define PAIR_STYLE

#include "field_component.h"

class PairStyle {
  private:
    int size ;
  public:
    int type1, type2 ;
    double *u, **f, **vir, energy, *total_vir ;
    double *rho1, *rho2, **force1, **force2 ;
    complex<double> *u_k, **f_k, **vir_k ;

    PairStyle( int, FieldComponent, FieldComponent ) ;
    ~PairStyle() ;
    
    double calc_energy() ;
    void setup_virial() ;
    void calc_all() ;
} ;

#endif
