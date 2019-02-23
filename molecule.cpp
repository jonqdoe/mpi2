#include "globals.h"
#include "molecule.h"

Molecule::Molecule( int sites_per_molec ) {

  this->Ns = sites_per_molec ;
  comp = ( int* ) calloc( sites_per_molec, sizeof(int) ) ;

}


Molecule::~Molecule() {

}
