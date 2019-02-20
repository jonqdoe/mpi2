#include "globals.h"


void ms_forces() {

  int i,j,k,m,n,o, id, gind, mod ;
  double ddot, factor = mu * gvol / rho0 ;
  
  double **T_tens ;
  T_tens = new double* [Dim] ;
  for ( m=0 ; m<Dim ; m++ ) {
    T_tens[m] = new double [Dim] ;
    for ( n=0 ; n<Dim ; n++ )
      T_tens[m][n] = 0.0 ;
  }

  // The first contribution to the force from
  // the M-S interactions. Arises from derivative
  // of \delta( r - r_i )
  for ( j=0 ; j<Dim ; j++ ) {

    // This loop builds up S(m,n)*\grad_j uG
    for ( m=0 ; m<Dim ; m++ ) {
      for ( n=m ; n<Dim ; n++ ) {
        
        // Take S[m][n] to k-space
        for ( i=0 ; i<ML ; i++ )
          tmp[i] = S_field[i][m][n] ;

        fftw_fwd( tmp, ktmp ) ;

        // Multiply by the gradient of uG in the j-direction:
        for ( i=0 ; i<ML ; i++ ) 
          ktmp[i] *= grad_uG_hat[j][i] ;

        // Back to real-space
        fftw_back( ktmp, tmp ) ;

        // Build up the tensor S(m,n) * \grad_j u_G
        for ( i=0 ; i<ML ; i++ ) {
          tmp_tensor[i][m][n] = tmp[i] ;
          if ( m != n )
            tmp_tensor[i][n][m] = tmp[i] ;
        }

      }
    }// for m=0:Dim-1



    // Apply the forces to particle i
    for ( i=0 ; i<ns_loc + total_ghost ; i++ ) {
      if ( i < ns_loc )
        id = my_inds[i] ;
      else
        id = ghost_inds[i-ns_loc] ;

      // Skip particles that don't carry orientation
      if ( !particle_has_orientation(id) )
        continue ;



      // Assign the forces to the particle
      // factor = mu * gvol / rho0
      for ( m=0 ; m<grid_per_partic ; m++ ) {

        gind = grid_inds[id][m] ;

        if ( gind == -1 )
          continue ;

        ddot = double_dot( mono_S[id], tmp_tensor[gind] ) ;

        f[id][j] -= ddot * grid_W[id][m] * factor ;

      }// for m=0:grid_per_partic-1

    }// for i=0:ns_loc+total_ghost-1
  }//for j=0:Dim-1


  // Second contribution to M-S force.
  // Arises from derivative of orientation u
  // w.r.t. r_j. 

  for ( j=0 ; j<Dim ; j++ ) {

    int j1, j2, j3 ;
    double dr21[Dim], dr32[Dim], mdr21, mdr21_sq, mdr32, mdr32_sq, dudrj[Dim] ;

    for ( i=0 ; i<ns_loc + total_ghost ; i++ ) {
      if ( i < ns_loc )
        id = my_inds[i] ;
      else
        id = ghost_inds[i-ns_loc] ;


      // Skip particles that don't carry orientation
      if ( !particle_has_orientation(id) )
        continue ;


      mod = id % ( Nda + Ndb ) ;
      j2 = id ;
      for ( m=0 ; m<Dim ; m++ ) 
        dudrj[m] = 0.0 ;

      ////////////////////////////////////////////////////////
      // First, build up vector du/dr_j2 in the j direction //
      ////////////////////////////////////////////////////////

      // Contribution from the particle before j2
      if ( mod != Nda ) {
        j1 = j2 - 1 ;
        mdr21_sq = pbc_mdr2( x[j2], x[j1], dr21 ) ;
        mdr21 = sqrt( mdr21_sq );
        
        for ( k=0 ; k<Dim ; k++ ) 
          dudrj[k] += 0.5 * ( KDelta(k,j) / mdr21 - dr21[k] * dr21[j] / mdr21  / mdr21_sq ) ;
      }

      if ( mod < ( Nda + Ndb -1 ) ) {
        j3 = j2 + 1 ;
        mdr32_sq = pbc_mdr2( x[j3], x[j2], dr32 ) ;
        mdr32 = sqrt(mdr32_sq) ;

        for ( m=0 ; m<Dim ; m++ ) 
          dudrj[m] += 0.5 * ( dr32[m] * dr32[j] / mdr32 / mdr32_sq - KDelta(m,j) / mdr32 ) ;

      }

      // Kill the pre-factor of 0.5 if this is first or last B monomer
      if ( mod == Nda || mod == ( Nda + Ndb - 1 ) ) {
        for ( m=0 ; m<Dim ; m++ )
          dudrj[m] *= 2.0 ;
      }


      //////////////////////////////////////
      // Form the tensor (Iu + uI)*dudr_j //
      //////////////////////////////////////
      for ( k=0 ; k<Dim ; k++ )
        for ( m=0 ; m<Dim ; m++ )
          T_tens[k][m] = 0.0 ;

      for ( k=0 ; k<Dim ; k++ ) 
        for ( m=0 ; m<Dim ; m++ ) 
          for ( n=0 ; n<Dim ; n++ ) 
            T_tens[m][n] += ( KDelta(k,m) * mono_u[j2][n] + mono_u[j2][k] * KDelta(m,n) ) * dudrj[n] ;


      //////////////////////////////
      // Dot T_tens with S_conv_u //
      //////////////////////////////
      
      // Assign the forces to the particle
      // factor = mu * gvol / rho0
      for ( m=0 ; m<grid_per_partic ; m++ ) {

        gind = grid_inds[id][m] ;

        if ( gind == -1 )
          continue ;

        ddot = double_dot( T_tens, S_conv_u[gind] ) ;

        f[id][j] -= ddot * grid_W[id][m] * factor ;

      }// for m=0:grid_per_partic-1
      

    }// for i=0:ns_local+total_ghost

  }// for j=0:Dim-1

  for ( n=0 ; n<Dim ; n++ )
    delete T_tens[n] ;
  delete T_tens ;

}





