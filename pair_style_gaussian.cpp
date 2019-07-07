#include "globals.h"
#include "pair_style_gaussian.h"
#include "field_component.h"


void Gaussian::Initialize_Gaussian(double Ao, double sigma_squared, int alloc_size, FieldComponent A, FieldComponent B) {
  Initialize_PairStyle(alloc_size, A, B) ;

  printf("Setting up Gaussian pair style...") ;fflush(stdout) ;
  int i , j;
  double k2, kv[Dim] ;


  // Define the potential and the force in k-space
  for ( i=0 ; i<alloc_size ; i++ ) {
    k2 = get_k(i, kv ) ;

    this->u_k[i] = Ao * exp( -k2 * sigma_squared / 2.0 ) ;
    for ( j=0 ; j<Dim ; j++ )
      this->f_k[j][i] = -I * kv[j] * this->u_k[i] ;

  }

  fftw_back( this->u_k, this->u ) ;
  for ( j=0 ; j<Dim ; j++ ) 
    fftw_back( this->f_k[j], this->f[j] ) ;

  this->setup_virial() ;

  printf("done!\n"); fflush(stdout) ;
}


void allocate_gaussians() {
  gaussian_prefactor = ( double* ) calloc( n_gaussian_pairstyles, sizeof(double) ) ;
  gaussian_sigma = ( double* ) calloc( n_gaussian_pairstyles, sizeof(double) ) ;
  gaussian_types = ( int** ) calloc( n_gaussian_pairstyles, sizeof( int* ) ) ;
  for (int i=0 ; i<n_gaussian_pairstyles ; i++ ) 
    gaussian_types[i] = (int*) calloc( 2, sizeof(int) ) ;
}


Gaussian::Gaussian() {

}

Gaussian::~Gaussian() {

}
