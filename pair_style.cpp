#include "globals.h"
#include "pair_style.h"


// Constructor, allocates memory for:
// potential, forces, virial contribution
// in both r- and k-space
void PairStyle::Initialize_PairStyle(int alloc_size, FieldComponent A, FieldComponent B ) {
  int n_off_diag = Dim + (Dim*Dim-Dim)/2 ;

  printf("setting up pair style!\n"); fflush(stdout) ;
  size = alloc_size ;

  this->u = ( double* ) calloc( alloc_size, sizeof(double) ) ;
  this->f = ( double** ) calloc( Dim, sizeof(double*) ) ;
  this->vir = ( double** ) calloc( n_off_diag, sizeof(double*) ) ;

  this->u_k = ( complex<double>* ) calloc(alloc_size, sizeof( complex<double> ) ) ;
  this->f_k = ( complex<double>** ) calloc(Dim, sizeof( complex<double>* ) ) ;
  this->vir_k = ( complex<double>** ) calloc(n_off_diag, sizeof( complex<double>* ) ) ;

  for ( int i=0 ; i<Dim ; i++ ) {
    this->f[i] = ( double* ) calloc( alloc_size, sizeof(double) ) ;
    this->f_k[i] = ( complex<double>* ) calloc( alloc_size, sizeof(complex<double>) ) ;
  }

  total_vir = ( double* ) calloc( n_off_diag, sizeof(double) ) ;

  // Vir stores the diagonal plus the off-diagonal terms
  // The term in parenthesis (Dim*Dim-Dim) will always be even
  for ( int i=0 ; i<n_off_diag ; i++ ) {
    this->vir[i] = ( double* ) calloc( alloc_size , sizeof(double) ) ;
    this->vir_k[i] = ( complex<double>* ) calloc( alloc_size , sizeof(complex<double>) ) ;
  }


  // Pointers to the two density fields involved in this interaction
  this->rho1 = A.rho ;
  this->rho2 = B.rho ;
  this->force1 = A.force ;
  this->force2 = B.force ;
  
  printf("Finished setting up pair style!\n"); fflush(stdout) ;
}




// Calculates the energy, virial, and gradU
// fields for this pair style.
void PairStyle::CalcAll() {
  
  int i, j ;
  fftw_fwd( this->rho2, ktmp ) ;

  // Energy calculation 
  for ( i=0 ; i<ML ; i++ ) 
    ktmp2[i] = ktmp[i] * this->u_k[i] ;

  fftw_back( ktmp2, tmp ) ;

  for ( i=0 ; i<ML ; i++ ) 
    tmp[i] *= this->rho1[i] ;

  this->energy = integrate(tmp) ;


  // Force calculation
  // rho2 acting on rho1
  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = ktmp[i] * this->f_k[j][i] ;

    fftw_back( ktmp2, tmp ) ;

    for ( i=0 ; i<ML ; i++ ) {
      this->force1[j][i] += this->rho1[i] * tmp[i] ;
    }
  }


  // Force calculation
  // rho1 acting on rho2
  fftw_fwd( this->rho1, ktmp ) ;
  for ( j=0 ; j<Dim ; j++ )  {

    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = ktmp[i] * this->f_k[j][i] ;
    
    fftw_back( ktmp2, tmp ) ;
    
    for ( i=0 ; i<ML ; i++ )
      this->force2[j][i] += this->rho2[i] * tmp[i] ;
  }

}






// Calculates the energy involved in this potential as
// energy = \int dr rho1(r) \int dr' u(r-r') rho2(r')
// The convolution theorem is used to efficiently evaluate 
// the integral over r'
double PairStyle::CalcEnergy( ) {
  
  int i ;
  fftw_fwd( this->rho2, ktmp ) ;

  for ( i=0 ; i<ML ; i++ ) 
    ktmp[i] *= this->u_k[i] ;

  fftw_back( ktmp, tmp ) ;

  for ( i=0 ; i<ML ; i++ ) 
    tmp[i] *= this->rho1[i] ;

  this->energy = integrate(tmp) ;
  return this->energy ;

}





// After an instance of a pair_style is initialized and
// u(r) and \mb f(r) is defined, this routine should be 
// called to set up the virial for this potential. The 
// form is general once the forces are defined.
void PairStyle::setup_virial() {
  int i, j ;
  double ro[Dim], dr[Dim], rv[Dim], mdr2 ;
  for ( j=0 ; j<Dim ; j++ ) 
    ro[j] = 0.0 ;

  // Calculate the virial function
  for ( i=0 ; i<this->size; i++ ) {
    get_r(i, rv) ;
    mdr2 = pbc_mdr2( rv, ro, dr ) ;
 
    for ( j=0 ; j<Dim ; j++ ) 
      this->vir[j][i] = dr[j] * this->f[j][i] ;
    
    // X-Y shear component
    this->vir[Dim][i] = dr[0] * this->f[1][i] ;

    // Add x-z, y-z shear components if relevant
    if ( Dim == 3 ) {
      this->vir[Dim+1][i] = dr[0] * this->f[2][i] ;
      this->vir[Dim+2][i] = dr[1] * this->f[2][i] ;
    }
  }
}





PairStyle::~PairStyle() {
  printf("PS here for some reason!\n"); fflush(stdout) ;
  free(this->u);
  free(this->u_k) ;

  for ( int i=0 ; i<ML ; i++ ) {
    free(this->f[i]) ;
    free(this->f_k[i]) ;
    free(this->vir[i]) ;
    free(this->vir_k[i]) ;
  }
  free(this->f) ;
  free(this->f_k) ;
  free(this->vir) ;
  free(this->vir_k) ;
}

PairStyle::PairStyle() {

}
