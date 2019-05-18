#include "globals.h"
#include "field_component.h"

void FieldComponent::Initialize( int alloc_size ) {

  this->rho = ( double* ) calloc( alloc_size, sizeof(double) ) ;
  this->force = ( double** ) calloc( Dim, sizeof( double* ) ) ;
  
  for ( int j=0 ; j<Dim ; j++ ) {
    this->force[j] = ( double* ) calloc( alloc_size, sizeof(double) ) ;
  }
  printf("entering zero_gradient..."); fflush(stdout);
  this->ZeroGradient() ;

  printf("DONE\n"); fflush(stdout);
}



void FieldComponent::ZeroGradient() {

  int i, j ;
  for ( j=0 ; j<Dim ; j++ )
    for ( i=0 ; i<ML ; i++ ) 
      this->force[j][i] = 0.0 ;
}



FieldComponent::FieldComponent( ) {
}

FieldComponent::~FieldComponent() {
}
