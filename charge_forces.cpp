#include "globals.h"

void charge_forces() {
 
  die("Charge forces not yet set up for MPI jobs!") ;

  //////////////////////////////
  // Solve Poisson's equation //
  //////////////////////////////

  fftw_fwd( rhoq, ktmp ) ; // F.T. the charge density
  for (int i=0 ; i<ML ; i++ ) {
    double k2, kv[Dim] ;

    k2 = get_k(i, kv) ;

    if ( k2 != 0.0 )
      ktmp2[i] = ktmp[i] * qhhat[i] * Lb / k2 ; //qhhat smears charge density

    else
      ktmp2[i] = 0.0 ;
  }//omp parallel for



  //////////////////////////////////
  // Save electrostatic potential //
  //////////////////////////////////
  fftw_back( ktmp2, psi ) ;


  //////////////////////////////
  // Calculate electric field //
  //////////////////////////////
  for ( int j=0 ; j<Dim ; j++ ) {

    for ( int i=0; i<ML ; i++ ) {
      double k2, kv[Dim] ;
      k2 = get_k( i , kv ) ;

      ktmp[i] = -I * kv[j] * ktmp2[i] ;
    }//end omp

    fftw_back( ktmp, Efield[j] ) ;
  }


  ///////////////////////////////////////////
  // Allocate Efield force on B components //
  ///////////////////////////////////////////

  // Only first B monomer charged //
  if ( Zb == 0.0 ) {
    for ( int i=Nda ; i<nsD ; i += (Nda+Ndb) ) {
      for ( int m=0 ; m<grid_per_partic ; m++ ) {

        int gind = grid_inds[i][m] ;

        for ( int j=0 ; j<Dim ; j++ ) 
          f[i][j] += Zb1 * Efield[j][ gind ] * grid_W[i][m] * gvol ;
      }//for (m=0...
    }//for ( i=Nda...
  }

  // All B monomers charged //
  else {
    for ( int i=0 ; i<nsD ; i++ ) {
      if ( tp[i] == 0 )
        continue ;
      
      for ( int m=0 ; m<grid_per_partic ; m++ ) {

        int gind = grid_inds[i][m] ;

        for ( int j=0 ; j<Dim ; j++ ) 
          f[i][j] += Zb * Efield[j][ gind ] * grid_W[i][m] * gvol ;
      }//for (m=0...

    }
  }

  // All C ions //
  for ( int i=nsD+nsA+nsB ; i<nsD+nsA+nsB+nsC ; i++ ) {
    if ( tp[i] == 0 )
      continue ;
    
    for ( int m=0 ; m<grid_per_partic ; m++ ) {

      int gind = grid_inds[i][m] ;

      for ( int j=0 ; j<Dim ; j++ ) 
        f[i][j] += Zc * Efield[j][ gind ] * grid_W[i][m] * gvol ;
    }//for (m=0...

  }


//  write_grid_data( "rhoq.dat", rhoq) ;
//  write_grid_data( "psi.dat", psi ) ;
//  write_grid_data( "Ex.dat", Efield[0] ) ;
//  write_grid_data( "Ey.dat", Efield[1] ) ;
//
//  die("writen!\n");
}
