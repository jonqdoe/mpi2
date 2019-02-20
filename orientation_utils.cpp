#include "globals.h"




// Takes the particle-level S tensors and adds them to the grid //
// It assumes that particle_orientations has been calculated    //
// prior to calling this function.
void add_S_to_grid( ) { 

  int i,j,k,m, id, Mid, mod ;
  double W3 ;

  
  for ( i=0 ; i<ns_loc+total_ghost ; i++ ) {
    if ( i < ns_loc ) 
      id = my_inds[i] ;
    else
      id = ghost_inds[ i-ns_loc ] ;

    if ( !particle_has_orientation( id ) )
      continue ;


    for ( m=0 ; m<grid_per_partic ; m++ ) {
      Mid = grid_inds[id][m] ;

      if ( Mid == -1 )
        continue ;

      W3 = grid_W[id][m] ;

      for ( j=0 ; j<Dim ; j++ )
        for ( k=0 ; k<Dim ; k++ ) 
          S_field[Mid][j][k] += W3 * mono_S[id][j][k] ;
    }// for m < grid_per_partic

  }// for i < ns_loc + total_ghost
}



// Calculates the u vector 
// for each diblock B monomer
void particle_orientations() {

  int i,j,id,k,m,n, mod,mod1, id1 ;
  double mdr, mdr2, dr[Dim] ;

  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;
    mod = id % ( Nda + Ndb ) ;

    if ( (id < nD * ( Nda + Ndb ))  // id is a diblock id
        && (mod >= Nda ) )           // id is a B monomer
    {
      
      // Reset monomer information
      for ( m=0 ; m<Dim ; m++ ) 
        mono_u[id][m] = 0.0 ;

      // Skip this part if id is a chain end
      if ( mod < (Nda + Ndb - 1) ) {
        id1 = id + 1 ;
 
        // Define unit vector
        mdr2 = pbc_mdr2( x[id], x[id1], dr ) ;
        mdr = sqrt( mdr2 ) ;
        for ( m=0 ; m<Dim ; m++ ) 
          dr[m] /= mdr ;
 
        // Accumulate u, S
        // Factor of 0.5 from averaging over involved bonds. Most monomers 
        // will participate in two bonds.
        for ( m=0 ; m<Dim ; m++ ) 
          mono_u[id][m]  += 0.5 * dr[m] ;
      }// if not a chain end


      // Skip this part if id is first B monomer
      if ( mod > Nda ) {
        id1 = id - 1 ;
 
        // Define unit vector
        mdr2 = pbc_mdr2( x[id1], x[id], dr ) ;
        mdr = sqrt( mdr2 ) ;
        for ( m=0 ; m<Dim ; m++ ) 
          dr[m] /= mdr ;
 
        // Accumulate u, S
        // Factor of 0.5 from averaging over involved bonds. Most monomers 
        // will participate in two bonds.
        for ( m=0 ; m<Dim ; m++ ) 
          mono_u[id][m]  += 0.5 * dr[m] ;
      }// if not first B monomer



      // Scale the magnitude of the terms for the first
      // B monomer and the chain end.
      if ( mod == Nda || mod == (Nda + Ndb - 1) ) {
        for ( m=0 ; m<Dim ; m++ ) 
          mono_u[id][m] *= 2.0 ;
      }


    }//if ( id is B part of diblock
  }// for i < ns_loc
}


// Calculates the S tensor for each particle.
// ASsumes u has been calculated and communicated
// in parallel across adjacent procs.
void particle_Stensor() {

  int i,j,id,k,m,n, mod,mod1, id1 ;
  double mdr, mdr2, dr[Dim] ;

  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;

    if ( !particle_has_orientation( id ) )
      continue ;


    for ( m=0 ; m<Dim ; m++ )
      for ( n=0 ; n<Dim ; n++ )
        mono_S[id][m][n] = mono_u[id][m] * mono_u[id][n] - KDelta(m,n)/Dim ;

  }// i=0:ns_loc
}



// Determines if a particle carries a u vector and
// and S tensor.
int particle_has_orientation( int id ) {

  int mod ;
  if ( id >= nD * ( Nda + Ndb ) )
    return 0 ;
  else if ( id % ( Nda + Ndb ) < Nda )
    return 0 ;
  else 
    return 1;
}


void make_eigenval_map() {

  int i,j ;
  double *Vec ;
  Vec = new double [Dim] ;

  for ( i=0 ; i<ML ; i++ ) {
    tmp[i] = diag_mat( S_field[i], Vec ) ;
  }

  delete Vec ;

  lc_order_param = integrate(tmp) / V ;

  write_grid_data( "eigen_val", tmp ) ;

}



// Calculates the convolution of the S field with the Gaussian potentail.
// This is needed for both force and energy evalulations. Should be 
// calc'd at the beginning of the force routine.
void calc_S_conv_uG() {

  int i,j,k,m,n ;

  fftw_fwd( uG, ktmp2 ) ;

  for ( m=0 ; m<Dim ; m++ ) 
    for ( n=0 ; n<Dim ; n++ ) {
      for ( i=0 ; i<ML ; i++ )
        tmp[i] = S_field[i][m][n] ;

      fftw_fwd( tmp, ktmp ) ;

      for ( i=0 ; i<ML ; i++ )
        ktmp[i] *= ktmp2[i] ;

      fftw_back( ktmp, tmp ) ;

      for ( i=0 ; i<ML ; i++ ) 
        S_conv_u[i][m][n] = tmp[i] ;

    }

}


