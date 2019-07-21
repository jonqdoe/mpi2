#include "globals.h"
#include "pair_style_erfc2.h"
#include "field_component.h"


void ErfcErfc::Initialize_Erfc2(double Ao, double Rp, double xi,
    int alloc_size, FieldComponent A, FieldComponent B) {

  Initialize_PairStyle(alloc_size, A, B) ;

  printf("Setting up Erfc-Erfc pair style...") ;fflush(stdout) ;
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
    ktmp[i] *= ktmp[i] ;

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

void allocate_erfc2s() {
  erfc2_prefactor = ( double* ) calloc( n_erfc2_pairstyles, sizeof(double) ) ;
  erfc2_Rp = ( double* ) calloc( n_erfc2_pairstyles, sizeof(double) ) ;
  erfc2_xi = ( double* ) calloc( n_erfc2_pairstyles, sizeof(double) ) ;
  erfc2_types = ( int** ) calloc( n_erfc2_pairstyles, sizeof( int* ) ) ;
  for (int i=0 ; i<n_erfc2_pairstyles ; i++ ) 
    erfc2_types[i] = (int*) calloc( 2, sizeof(int) ) ;
}


ErfcErfc::ErfcErfc() {

}

ErfcErfc::~ErfcErfc() {

}


