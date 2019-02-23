#include "globals.h"
#include "pair_style_gaussian.h"
#include "field_component.h"


Gaussian::Gaussian(double Ao, double sigma, int alloc_size, FieldComponent A, FieldComponent B) :
  PairStyle(alloc_size, A, B) {

  int i , j;
  double k2, kv[Dim] ;


  // Define the potential and the force in k-space
  for ( i=0 ; i<alloc_size ; i++ ) {
    k2 = get_k(i, kv ) ;

    this->u_k[i] = Ao * exp( -k2 * sigma / 2.0 ) ;
    for ( j=0 ; j<Dim ; j++ )
      this->f_k[j][i] = I * kv[j] * this->u_k[i] ;

  }

  fftw_back( this->u_k, this->u ) ;
  for ( j=0 ; j<Dim ; j++ ) 
    fftw_back( this->f_k[j], this->f[j] ) ;

  this->setup_virial() ;

}

