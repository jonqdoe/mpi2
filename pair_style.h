#ifndef PAIR_STYLE
#define PAIR_STYLE
#include <complex>
#include "field_component.h"
using namespace std ;

class PairStyle {
  private:
    int size ;
  public:
    int type1, type2 ;
    double *u, **f, **vir, energy, *total_vir ;
    double *rho1, *rho2, **force1, **force2 ;
    complex<double> *u_k, **f_k, **vir_k ;

    PairStyle( ) ;
    ~PairStyle() ;
    
    void Initialize_PairStyle( int, FieldComponent, FieldComponent ) ;
    double CalcEnergy() ;
    void setup_virial() ;
    void CalcAll() ;
} ;

#endif
