#include "globals.h"
#include "pair_style_gauss-erfc.h"
#include "field_component.h"


void GaussianErfc::Initialize_GaussianErfc(double Ao, double sigma_squared, 
    double Rp, double xi, int alloc_size, FieldComponent A, FieldComponent B) {

  Initialize_PairStyle(alloc_size, A, B) ;

  printf("Setting up Gaussian-Erfc pair style...") ;fflush(stdout) ;
  int i , j;
  double k2, kv[Dim], ro[Dim], rc[Dim], dr[Dim], mdr2, mdr ;

  for ( j=0 ; j<Dim ; j++ )
    ro[j] = 0.0 ;


  // Define the potential and the force in k-space
  for ( i=0 ; i<alloc_size ; i++ ) {
    get_r( i, rc ) ;

    mdr2 = pbc_mdr2( ro, rc, dr) ;
    mdr = sqrt(mdr2) ;

    tmp[i] = Ao * V / 2.0 * ( 1.0 - erf( ( mdr - Rp ) / xi ) ) ;
  }

  fftw_fwd( tmp, ktmp ) ;


  for ( i=0 ; i<ML ; i++ ) {
    k2 = get_k(i, kv) ;

    ktmp[i] *= exp( -k2 * sigma_squared / 2.0 ) ;

    this->u_k[i] = ktmp[i] ;

    for ( j=0 ; j<Dim ; j++ )
      this->f_k[j][i] = -I * kv[j] * this->u_k[i] ;
  }


  fftw_back( this->u_k, this->u ) ;
  for ( j=0 ; j<Dim ; j++ ) 
    fftw_back( this->f_k[j], this->f[j] ) ;

  this->setup_virial() ;

  printf("done!\n"); fflush(stdout) ;
}

void allocate_gausserfcs() {
  gausserfc_prefactor = ( double* ) calloc( n_gausserfc_pairstyles, sizeof(double) ) ;
  gausserfc_sigma = ( double* ) calloc( n_gausserfc_pairstyles, sizeof(double) ) ;
  gausserfc_Rp = ( double* ) calloc( n_gausserfc_pairstyles, sizeof(double) ) ;
  gausserfc_xi = ( double* ) calloc( n_gausserfc_pairstyles, sizeof(double) ) ;
  gausserfc_types = ( int** ) calloc( n_gausserfc_pairstyles, sizeof( int* ) ) ;
  for (int i=0 ; i<n_gausserfc_pairstyles ; i++ ) 
    gausserfc_types[i] = (int*) calloc( 2, sizeof(int) ) ;
}


GaussianErfc::GaussianErfc() {

}

GaussianErfc::~GaussianErfc() {

}


