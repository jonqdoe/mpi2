#include "globals.h"
#include "timing.h"

void spline_get_weights( double , double , double* ) ;
void add_segment( int ) ;
void zero_fields( void ) ;



///////////////////////////////////
// Add all particles to the grid //
///////////////////////////////////
void charge_grid( ) {

  grid_t_in = time(0) ;

  int i, j, k ;

  zero_fields( ) ;

  ////////////////////////////////////////////////////
  // Add segments to each processors density fields //
  ////////////////////////////////////////////////////
  // Add full-bodied particles //
  time_debug_in = time(0);
  for ( i=0 ; i<ns_loc ; i++ ) {
    int id = my_inds[i] ;
    add_segment( id ) ;
  }
  time_debug_out = time(0);
  time_debug_tot_time += time_debug_out - time_debug_in ;

  // Add ghost particles
  for ( i=0 ; i<total_ghost ; i++ ) {
    int id = ghost_inds[i] ;
    add_segment( id ) ;
  }


  grid_t_out = time(0) ;
  grid_tot_time += ( grid_t_out - grid_t_in ) ;
}




//////////////////////////////////////////////
// Adds the segment associated with particle //
// "id" to the PME grid using Lagrange      //
// interpolation scheme. From JCP V103 3668 //
//////////////////////////////////////////////
void add_segment( int id ) {



  int j, g_ind[Dim] , ix, iy, iz, nn[Dim] , Mindex, grid_ct; 
  double **W , gdx , W3;
  
  W = ( double** ) calloc( Dim , sizeof( double* ) );



  ///////////////////////////////////////////////
  // First, determine the relevant weights for //
  // all grid points in all directions.        //
  ///////////////////////////////////////////////
  for ( j=0 ; j<Dim ; j++ ) {

    W[j] = ( double* ) calloc( pmeorder+1 , sizeof( double ) );


    // Distance to nearest grid point if even //
    if ( pmeorder % 2 == 0 ) {
      g_ind[j] = int( ( x[id][j] + 0.5 * dx[j] ) / dx[j] ) ;

      gdx = x[id][j] - double( g_ind[j] ) * dx[j] ;
    }
 

    // Distance to nearest mid-point between grid points if odd //
    else {
      g_ind[j] = int( ( x[id][j]  ) / dx[j] ) ;

      if ( g_ind[j] >= Nx[j] ) {
        char last_words[80] ;
        sprintf( last_words, "Particle %d out of bounds: %lf %lf %lf\n",
            id, x[id][0], x[id][1], x[id][2] ) ;
        die(last_words) ;
      }

      gdx = x[id][j]  - ( double( g_ind[j] ) + 0.5 ) * dx[j] ;
    }


    /////////////////////////////////////////
    // Get the weights for each grid point //
    /////////////////////////////////////////
    spline_get_weights( gdx , dx[j] , W[j] );

  }//for ( j=0 ; j<3...




  ////////////////////////////////////////////////////
  // Assign the weights to all relevant grid points //
  ////////////////////////////////////////////////////
  grid_ct = 0 ;
  
  ///////////////////////////////////////////
  // 3D version of particle-to-mesh scheme //
  ///////////////////////////////////////////
  if ( Dim == 3 ) {
    int order_shift = pmeorder/2 ;

    for ( ix = 0 ; ix < pmeorder+1 ; ix++ ) {
      
      nn[0] = g_ind[0] + ix - order_shift ;
      
      if ( nn[0] < 0 ) nn[0] += Nx[0] ;
      else if ( nn[0] >= Nx[0] ) nn[0] -= Nx[0] ;
  
      for ( iy = 0 ; iy < pmeorder+1 ; iy++ ) {
  
        nn[1] = g_ind[1] + iy - order_shift ;
        
        if ( nn[1] < 0 ) nn[1] += Nx[1] ;
        else if ( nn[1] >= Nx[1] ) nn[1] -= Nx[1] ;
  
        for ( iz = 0 ; iz < pmeorder+1 ; iz++ ) {
  
          nn[2] = g_ind[2] + iz - order_shift ;
          
          if ( nn[2] < 0 ) nn[2] += Nx[2] ;
          else if ( nn[2] >= Nx[2] ) nn[2] -= Nx[2] ;
  

          // Ensure this z-grid is on the current processor
          if ( nn[2] < zstart || nn[2] >= zstart + NxL[2] ) {
            grid_inds[ id ][ grid_ct ] =  -1 ;
            grid_W[ id ][ grid_ct ] = 0.0 ;
  
            grid_ct++ ;

            continue ;
          }
 

          // stack_to_local() returns index in [0,ML]
          Mindex = stack_to_local( nn ) ;
  
          if ( Mindex >= ML ) {
            char kill_msg[80];
            sprintf( kill_msg, "Mindex = %d out of range! nn = [%d, %d]\n", Mindex, nn[0], nn[1] ) ;
            die(kill_msg);
          }
  
          W3 = W[0][ix] * W[1][iy] * W[2][iz] / gvol ;
 

          Components[ tp[id] ].rho[Mindex] += W3 ;


          grid_inds[ id ][ grid_ct ] = Mindex ;
          grid_W[ id ][ grid_ct ] = W3 ;
  
          grid_ct++ ;
  
        }
      }
    }
  }

  ///////////////////////////////////////////
  // 2D version of particle-to-mesh scheme //
  ///////////////////////////////////////////
  else if ( Dim == 2 ) {
    int order_shift = pmeorder/2 ;

    for ( ix = 0 ; ix < pmeorder+1 ; ix++ ) {
      
      nn[0] = g_ind[0] + ix - order_shift ;
      
      if ( nn[0] < 0 ) nn[0] += Nx[0] ;
      else if ( nn[0] >= Nx[0] ) nn[0] -= Nx[0] ;
  
      for ( iy = 0 ; iy < pmeorder+1 ; iy++ ) {
  
        nn[1] = g_ind[1] + iy - order_shift ;
        
        if ( nn[1] < 0 ) nn[1] += Nx[1] ;
        else if ( nn[1] >= Nx[1] ) nn[1] -= Nx[1] ;
 


        // Ensure we're working with a grid location on the current processor
        if ( nn[1] < zstart || nn[1] >= zstart + NxL[1] ) {
          grid_inds[ id ][ grid_ct ] =  -1 ;
          grid_W[ id ][ grid_ct ] = 0.0 ;
        
          grid_ct++ ;
          continue ;
        }

        Mindex = stack_to_local( nn ) ;
  
        if ( Mindex >= ML ) {
          char kill_msg[80];
          sprintf( kill_msg, "Mindex = %d out of range! nn = [%d, %d] NxL = [%d, %d]\n", 
              Mindex, nn[0], nn[1], NxL[0], NxL[1] ) ;
          die(kill_msg);
        }

  
        W3 = W[0][ix] * W[1][iy] / gvol ;
  
        
        Components[ tp[id] ].rho[Mindex] += W3 ;

  
        grid_inds[ id ][ grid_ct ] = Mindex ;
        grid_W[ id ][ grid_ct ] = W3 ;
  
        grid_ct++ ;
  
        
      }
    }
  }


  ///////////////////////////
  // Free allocated memory //
  ///////////////////////////
  for ( j=0 ; j<Dim ; j++ ) 
    free( W[j] ) ;

  free(W) ;

}// End lagrange add charge








///////////////////////////////////////////////////////
// Taken from Deserno and Holm, JCP V109 7678 (1998) //
///////////////////////////////////////////////////////
void spline_get_weights( double dx , double H , double *W ) {
  double sx = dx / H ;
  
  double sx2, sx3, sx4, sx5;
  double scale = double( M ) / V ;

  if ( pmeorder == 0 ) 
    W[0] = 1. * scale ;

  else if ( pmeorder == 1 ) {
    W[0] = ( 0.5 - sx ) ;
    W[1] = ( 0.5 + sx ) ;
  }

  else if ( pmeorder == 2 ) {
    sx2 = sx * sx ;

    W[0] = (0.125 - 0.5 * sx + 0.5 * sx2)  ;
    W[1] = (0.75 - sx2 ) ;
    W[2] = (0.125 + 0.5 * sx + 0.5*sx2 ) ;

  }

  else if ( pmeorder == 3 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx ;

    W[0] = ( 1. - 6.*sx + 12.*sx2 - 8.*sx3 ) / 48. ;
    W[1] = ( 23. - 30.*sx - 12.*sx2 + 24.*sx3 ) / 48. ;
    W[2] = ( 23. + 30.*sx - 12.*sx2 - 24.*sx3 ) / 48. ;
    W[3] = ( 1. + 6.*sx + 12.*sx2 + 8.*sx3 ) / 48. ;
  }

  else if ( pmeorder == 4 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx2 ;
    sx4 = sx2 * sx2 ;

    W[0] = (1. - 8.*sx + 24.*sx2 - 32.*sx3 + 16.*sx4 ) / 384.;
    W[1] = (19. - 44.*sx + 24.*sx2 + 16.*sx3 - 16.*sx4 ) / 96. ;
    W[2] = (115. - 120.*sx2 + 48.*sx4 ) / 192. ;
    W[3] = (19. + 44.*sx + 24.*sx2 - 16.*sx3 - 16.*sx4 ) / 96. ;
    W[4] = (1. + 8.*sx + 24.*sx2 + 32.*sx3 + 16.*sx4 ) / 384.;
  }

  else if ( pmeorder == 5 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx2 ;
    sx4 = sx2 * sx2 ;
    sx5 = sx4 * sx ;

    W[0] = (1. - 10.*sx + 40.*sx2 - 80.*sx3 + 80.*sx4 - 32.*sx5 ) / 3840.;
    W[1] = (237. - 750.*sx + 840.*sx2 - 240.*sx3 - 240.*sx4 + 160.*sx5 ) / 3840.;
    W[2] = (841. - 770.*sx - 440.*sx2 + 560.*sx3 + 80.*sx4 - 160.*sx5 ) / 1920. ;
    W[3] = (841. + 770.*sx - 440.*sx2 - 560.*sx3 + 80.*sx4 + 160.*sx5 ) / 1920. ;
    W[4] = (237. + 750.*sx + 840.*sx2 + 240.*sx3 - 240.*sx4 - 160.*sx5 ) / 3840.;
    W[5] = (1. + 10.*sx + 40.*sx2 + 80.*sx3 + 80.*sx4 + 32.*sx5 ) / 3840.;

  }


  else
    die("P3M not set up for this interpolation order!\n");

}


void zero_fields() {
  int i,j,k ;
 
  for ( i=0 ; i<ntypes ; i++ ) {
    for ( j=0 ; j<ML ; j++ )
      Components[i].rho[j] = 0.0 ;
  }

}
