#include "globals.h"
#include "pair_style.h"


// Constructor, allocates memory for:
// potential, forces, virial contribution
// in both r- and k-space
PairStyle::PairStyle(int alloc_size, double *rhoA, double *rhoB) {
  u = ( double* ) calloc( alloc_size, sizeof(double) ) ;
  f = ( double** ) calloc( alloc_size, sizeof(double*) ) ;
  vir = ( double** ) calloc( alloc_size, sizeof(double*) ) ;

  u_k = ( complex<double>* ) calloc(alloc_size, sizeof( complex<double> ) ) ;
  f_k = ( complex<double>** ) calloc(alloc_size, sizeof( complex<double>* ) ) ;
  vir_k = ( complex<double>** ) calloc(alloc_size, sizeof( complex<double>* ) ) ;

  for ( int i=0 ; i<alloc_size ; i++ ) {
    f[i] = ( double* ) calloc( Dim, sizeof(double) ) ;
    f_k[i] = ( complex<double>* ) calloc( Dim, sizeof(complex<double>) ) ;

    // Vir stores the diagonal plus the off-diagonal terms
    // The term in parenthesis (Dim*Dim-Dim) will always be even
    vir[i] = ( double* ) calloc( Dim + (Dim*Dim-Dim)/2 , sizeof(double) ) ;
    vir_k[i] = ( complex<double>* ) calloc( Dim + (Dim*Dim-Dim)/2 , sizeof(complex<double>) ) ;
  }


  // Pointers to the two density fields involved in this interaction
  rho1 = rhoA ;
  rho2 = rhoB ;
}







// Calculates the energy involved in this potential as
// energy = \int dr rho1(r) \int dr' u(r-r') rho2(r')
// The convolution theorem is used to efficiently evaluate 
// the integral over r'
double PairStyle::calc_energy( ) {
  
  int i ;
  fftw_fwd( this->rho2, ktmp ) ;

  for ( i=0 ; i<ML ; i++ ) 
    ktmp[i] *= this->u_k[i] ;

  fftw_back( ktmp, tmp ) ;

  for ( i=0 ; i<ML ; i++ ) 
    tmp[i] *= this->rho1[i] ;

  return integrate( tmp ) ;

}




PairStyle::~PairStyle() {
  free(u);
  free(u_k) ;

  for ( int i=0 ; i<ML ; i++ ) {
    free(f[i]) ;
    free(f_k[i]) ;
    free(vir[i]) ;
    free(vir_k[i]) ;
  }
  free(f) ;
  free(f_k) ;
  free(vir) ;
  free(vir_k) ;
}
