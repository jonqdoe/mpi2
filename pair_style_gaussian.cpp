#include "globals.h"
#include "pair_style_gaussian.h"


Gaussian::Gaussian(double Ao, double sigma, int alloc_size, double *rho1, double *rho2):
  PairStyle(alloc_size, rho1, rho2) {

  int i ;
  double k2, kv[Dim] ;

  for ( i=0 ; i<alloc_size ; i++ ) {
    k2 = get_k(i, kv ) ;

    this->u_k[i] = Ao * exp( -k2 * sigma / 2.0 ) ;
  }

  fftw_back( this->u_k, this->u ) ;
}

