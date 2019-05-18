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
 

  for ( j=0 ; j<ntypes ; j++ ) 
    Components[j].ZeroGradient() ;


  
  /////////////////////////////////////
  // Calculate the grid-based forces //
  /////////////////////////////////////
  for ( j=0 ; j<n_gaussian_pairstyles ; j++ )
    Gauss[j].CalcAll() ;




      

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
        // This division by rho(r) might be better off pre-calculated
        f[id][j] += Components[ tp[id] ].force[j][gind] / Components[ tp[id] ].rho[gind] ;

        //if ( tp[id] == 0 )
        //  f[id][j] -= gradwA[ j ][ gind ] * grid_W[id][m] ;
        //else if ( tp[id] == 1 )
        //  f[id][j] -= gradwB[ j ][ gind ] * grid_W[id][m] ;
        //else if ( tp[id] == 2 )
        //  f[id][j] -= gradwP[ j ][ gind ] * grid_W[id][m] ;
        //else if ( tp[id] == 3 )
        //  f[id][j] -= gradwC[ j ][ gind ] * grid_W[id][m] ;
      }
    }

    for ( j=0 ; j<Dim ; j++ )
      f[id][j] *= gvol ;
  }
 

 
  cout << "Entering bonds..." ;
  ////////////////////////////
  // Call the bonded forces //
  ////////////////////////////
  if ( n_total_bonds > 0 )
    bonds() ;

  cout << "done!\n" ;

  if ( n_total_angles > 0 ) 
    angles() ;




  if ( nprocs > 1 )
    communicate_ghost_forces() ;

}
