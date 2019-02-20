#include "globals.h"

void angles( void ) {

  int i,j, k, i1, i2, i3 ;
  double rij[Dim], rkj[Dim], mdrij, mdrkj, cos_th, dot, u_loc,
         mdrij2, mdrkj2 ;

  Uangle = u_loc = 0.0 ;

//  return ;

  for ( i=0 ; i<ns_loc ; i++ ) {
    i1 = my_inds[i] ;

    for ( j=0 ; j<n_angles[i1] ; j++ ) {

      if ( i1 != angle_first[i1][j] )
        continue ;

      i2 = angle_mid[i1][j] ;
      i3 = angle_end[i1][j] ;

      //int find_sum = 0 ;
      //find_sum += local_flag[i2] + is_partic_in_list( total_ghost, i2, ghost_inds ) ;
      //find_sum += local_flag[i3] + is_partic_in_list( total_ghost, i3, ghost_inds ) ;
      //if ( find_sum < 2 ) printf("Particles i2=%d or i3=%d missing in angle with i1=%d",i2,i3,i1);

      // Generate the vectors rij, rkj
      mdrij2 = pbc_mdr2( x[i1] , x[i2] , rij );
      mdrkj2 = pbc_mdr2( x[i3] , x[i2] , rkj );
 
      // Take dot product and deteremine the magnitude of rij, rkj
      dot = 0.0 ;
 
      for ( k=0 ; k<Dim ; k++ ) 
        dot += rij[k] * rkj[k] ;
 
      mdrij = sqrt( mdrij2 );
      mdrkj = sqrt( mdrkj2 );
 
      // cosine of angle between vectors
      cos_th = dot / mdrij / mdrkj ;


      u_loc += angle_coeff[i1][j] * ( 1.0 + cos_th ) ;


      for ( k=0 ; k<Dim ; k++ ) {
        double fi = angle_coeff[i1][j] * 
                    ( cos_th * rij[k] / mdrij2 - rkj[k] / mdrij / mdrkj ) ;
        double fk = angle_coeff[i1][j] * 
                    ( cos_th * rkj[k] / mdrkj2 - rij[k] / mdrij / mdrkj ) ;
        f[i1][k] += fi ;
        f[i2][k] += -( fi + fk ) ;
        f[i3][k] += fk ;
      }


    }//for j < n_angles[i1]

  }//for ( i < ns_loc

  MPI_Allreduce( &u_loc, &Uangle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;

}




void angle_init( void ) {

  int i, m, ind ;

  for ( i=0 ; i<nstot ; i++ ) 
    n_angles[i] = 0 ;


  // Rod block bonds //
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=Nda ; m<Nda + Ndb - 1 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;
      n_bonds[ind] = 0 ;

      bonded_to[ ind ][ n_bonds[ind] ] = ind+1 ;

      bond_eq[ ind ][ n_bonds[ind] ] = semiflex_req ;
      bond_coeff[ ind ][ n_bonds[ind] ] = semiflex_k ;

      n_bonds[ind] = 1 ;
    }
  } // for ( i=0 ; i<nT[k]



  // Rod block angles //
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=Nda ; m<Nda + Ndb - 2 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      angle_first[ ind ][ n_angles[ind] ] = ind ;
      angle_mid[ ind ][ n_angles[ind] ] = ind+1 ;
      angle_end[ ind ][ n_angles[ind] ] = ind+2 ;
      angle_coeff[ ind ][ n_angles[ind] ] = semiflex_lam ;
      n_angles[ind]++ ;
      
      
      angle_first[ ind+1 ][ n_angles[ind+1] ] = ind ;
      angle_mid[ ind+1 ][ n_angles[ind+1] ] = ind+1 ;
      angle_end[ ind+1 ][ n_angles[ind+1] ] = ind+2 ;
      angle_coeff[ ind+1 ][ n_angles[ind+1] ] = semiflex_lam ;
      n_angles[ind+1]++ ;
      

      angle_first[ ind+2 ][ n_angles[ind+2] ] = ind ;
      angle_mid[ ind+2 ][ n_angles[ind+2] ] = ind+1 ;
      angle_end[ ind+2 ][ n_angles[ind+2] ] = ind+2 ;
      angle_coeff[ ind+2 ][ n_angles[ind+2] ] = semiflex_lam ;
      n_angles[ind+2]++ ;
    }
  } // for ( i=0 ; i<nT[k]

}
