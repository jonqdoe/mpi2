#include "globals.h"
#include "field_component.h"

FieldComponent::FieldComponent( int alloc_size ) {

  this->rho = ( double* ) calloc( alloc_size, sizeof(double) ) ;
  this->gradU = ( double** ) calloc( Dim, sizeof( double* ) ) ;
  
  for ( int j=0 ; j<Dim ; j++ ) {
    this->gradU[j] = ( double* ) calloc( alloc_size, sizeof(double) ) ;
  }
}



FieldComponent::~FieldComponent() {
//  printf("FC here for some reason!\n"); fflush(stdout) ;
//  free(rho) ;
//  for ( int j=0 ; j<Dim ; j++ )
//    free( this->gradU[j] ) ;
//
//  free(this->gradU) ;
}
