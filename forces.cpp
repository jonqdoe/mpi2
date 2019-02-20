#include "globals.h"
void charge_grid( void ) ;
void bonds( void ) ;
void charge_forces( void ) ;
void communicate_ghost_forces( void ) ;
void angles( void ) ;
void particle_orientations( void ) ;
void calc_S_conv_uG(void);
void ms_forces( void ) ;



void forces() {

  int i,j, m, gind, t1, t2, id ;


  charge_grid() ;


  ///////////////////////////////////////////////
  // Reset the particle forces and grid grad w //
  ///////////////////////////////////////////////
  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;
    for ( j=0 ; j<Dim ; j++ )
      f[id][j] = 0.0 ;
  }
  for ( i=0 ; i<total_ghost ; i++ ) {
    id = ghost_inds[i] ;
    for ( j=0 ; j<Dim ; j++ )
      f[id][j] = 0.0 ;
  }
 

  for ( i=0 ; i<ML ; i++ )
    for ( j=0 ; j<Dim ; j++ ) 
      gradwA[j][i] = gradwB[j][i] = gradwC[j][i] = gradwP[j][i] = 0.0 ;



  //////////////////////////////////////////////////
  // Accumulate the monomer-monomer contributions //
  //////////////////////////////////////////////////
  
  // A acting on B, C //
  fftw_fwd( rho[0] , ktmp ) ;

  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

    for ( i=0 ; i<ML ; i++ ) {
      if ( chiAB != 0.0 )
        gradwB[j][i] += tmp[i] * chiAB / rho0 ;
      if ( chiAC != 0.0 )
        gradwC[j][i] += tmp[i] * chiAC / rho0 ;
    }
  }


  // B acting on A, C //
  fftw_fwd( rho[1] , ktmp ) ;
  
  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

    for ( i=0 ; i<ML ; i++ ) {
      if ( chiAB != 0.0 )
        gradwA[j][i] += tmp[i] * chiAB / rho0 ;
      if ( chiBC != 0.0 )
        gradwC[j][i] += tmp[i] * chiBC / rho0 ;
    }
  }

  // C acting on A, B //
  fftw_fwd( rho[3] , ktmp ) ;
  
  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

    for ( i=0 ; i<ML ; i++ ) {
      if ( chiAC != 0.0 )
        gradwA[j][i] += tmp[i] * chiAC / rho0 ;
      if ( chiBC != 0.0 )
        gradwB[j][i] += tmp[i] * chiBC / rho0 ;
    }
  }




  // Compressibility contribution //
  for ( i=0 ; i<ML ; i++ )
    tmp[i] = rhot[i] ;

  fftw_fwd( tmp , ktmp ) ;

  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

    for ( i=0 ; i<ML ; i++ ) {
      gradwA[j][i] += tmp[i] * kappa / rho0 ;
      gradwB[j][i] += tmp[i] * kappa / rho0 ;
      gradwC[j][i] += tmp[i] * kappa / rho0 ;
    }

  }


  ///////////////////////////////////////////
  // Accumulate the particle contributions //
  ///////////////////////////////////////////
  if ( nP > 0 ) {

    // Particle-particle //
    fftw_fwd( rho[2] , ktmp ) ;

    for ( j=0 ; j<Dim ; j++ ) {
      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uP_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * kappa / rho0 ;
    }

    // Particles acting on monomers //
    for ( j=0 ; j<Dim ; j++ ) {
      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ ) {
        gradwA[j][i] += tmp[i] * ( kappa + ( B_partics ? chiAB : 0.0 ) ) / rho0 ;
        gradwB[j][i] += tmp[i] * kappa ;
        gradwC[j][i] += tmp[i] * ( kappa + ( B_partics ? chiBC : 0.0 ) ) / rho0 ;
      }
    }

    // Add attractive contribution //
    if ( eps != 0.0 ) {
      for ( j=0 ; j<Dim ; j++ ) {

        for ( i=0 ; i<ML ; i++ )
          ktmp2[i] = grad_uAG_hat[j][i] * ktmp[i] ;
        fftw_back( ktmp2, tmp ) ;
       
        for ( i=0 ; i<ML ; i++ ) 
          gradwA[j][i] -= tmp[i] * eps / rho0 ;

      }// j=0:Dim-1
    }// if ( eps!= 0.0)


    // A MLonomers acting on particles //
    fftw_fwd( rho[0] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++ ) {

      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( B_partics ? chiAB : 0.0 ) ) / rho0 ;
    }

    if ( eps != 0.0 ) {
      for ( j=0 ; j<Dim ; j++ ) {
        for ( i=0 ; i<ML ; i++ )
          ktmp2[i] = grad_uAG_hat[j][i] * ktmp[i] ;

        fftw_back( ktmp2, tmp ) ;

        for ( i=0 ; i<ML ; i++ )
          gradwP[j][i] -= tmp[i] * eps/rho0 ;

      } // j=0:Dim-1
    }// eps != 0.0

    // B MLonomers acting on particles //
    fftw_fwd( rho[1] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++) {

      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * kappa / rho0 ;
    }

    // C Monomers acting on particles //
    fftw_fwd( rho[3] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++) {

      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( B_partics ? chiBC : 0.0 ) ) / rho0 ;
    }


  }// if ( nP > 0 )



      

  //////////////////////////////////////////////////////
  // Accumulate the nonbonded forces on each particle //
  //////////////////////////////////////////////////////
  for ( i=0 ; i<ns_loc+total_ghost ; i++ ) {

    if ( i < ns_loc )
      id = my_inds[i] ;
    else
      id = ghost_inds[ i-ns_loc ] ;


    for ( m=0 ; m < grid_per_partic ; m++ ) {

      // Note that in charge_grid(), grid_inds is populated with 
      // indices in the range [0,ML]. 
      // if gind == -1, then grid location is skipped.
      gind = grid_inds[ id ][ m ] ;


      if ( gind == -1 )
        continue ;


      for ( j=0 ; j<Dim ; j++ ) {
        if ( tp[id] == 0 )
          f[id][j] -= gradwA[ j ][ gind ] * grid_W[id][m] ;
        else if ( tp[id] == 1 )
          f[id][j] -= gradwB[ j ][ gind ] * grid_W[id][m] ;
        else if ( tp[id] == 2 )
          f[id][j] -= gradwP[ j ][ gind ] * grid_W[id][m] ;
        else if ( tp[id] == 3 )
          f[id][j] -= gradwC[ j ][ gind ] * grid_W[id][m] ;
      }
    }

    for ( j=0 ; j<Dim ; j++ )
      f[id][j] *= gvol ;
  }

 

  ////////////////////////////
  // Call the charge forces // 
  ////////////////////////////
  if ( Lb > 0.0 )
    charge_forces() ;
  
 
  ////////////////////////////
  // Call the bonded forces //
  ////////////////////////////
  bonds() ;

  if ( semiflex )
    angles() ;


  
  if ( mu != 0.0 ) {
    calc_S_conv_uG() ;
    ms_forces() ;
  }


  if ( nprocs > 1 )
    communicate_ghost_forces() ;

}
