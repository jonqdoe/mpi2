#include "globals.h"
#include "molecule_homopolymer.h"

Homopolymer::Homopolymer(int sites_per_molec, int component ) 
  : Molecule( sites_per_molec ) {

    int i ;
    for ( i=0 ; i<sites_per_molec ; i++ ) 
      comp[i] = component ;

}
