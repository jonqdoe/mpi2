#include "globals.h"


void bonds() {

  int i,j, k, id1, id2, present_flag, btype ;
  double mdr2, mdr, dr[Dim], delr, mf ;
  double ub_loc = 0.0 , utmp ;


  for ( i=0 ; i<ns_loc ; i++ ) {
    id1 = my_inds[i] ;

    for ( j=0 ; j<n_bonds[id1] ; j++ ) {
      id2 = bonded_to[ id1 ][ j ] ;

      if ( id2 < id1 )
        continue ;

      btype = bond_type[id1][j] ;

      mdr2 = pbc_mdr2( x[id1], x[id2], dr ) ;

      mdr = sqrt(mdr2) ;

//      cout << j << " " << id2 << " " << btype << endl;

      delr = mdr - bond_eq[ btype ] ;

      ub_loc += delr * delr * bond_coeff[btype] / 2.0 ;
      
      mf = bond_coeff[btype] * delr / mdr ;
      
      for ( k=0 ; k<Dim ; k++ ) {
        f[id1][k] -= mf * dr[k] ;
        f[id2][k] += mf * dr[k] ;
      }

    }//for ( j=0 ; j<n_bonds ;
    die("through i=0");
  }// for ( i=0 ; i<ns_loc


  MPI_Allreduce( &ub_loc, &utmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;
  Ubond = utmp ;

}




void bond_stress() {

  int i,j, k1, k2, id1, id2 ;
  int btype ;
  double mdr2, mdr, dr[Dim], delr, mf ;

  for ( i=0 ; i<ns_loc ; i++ ) {
    id1 = my_inds[i] ;

    for ( j=0 ; j<n_bonds[id1] ; j++ ) {
      id2 = bonded_to[ id1 ][ j ] ;

      if ( id2 < id1 ) 
        continue ;

      btype = bond_type[id1][j] ;

      mdr2 = pbc_mdr2( x[id1], x[id2], dr ) ;

      mdr = sqrt(mdr2) ;

      delr = mdr - bond_eq[ btype ] ;
      
      mf = bond_coeff[btype] * delr / mdr ;
      for ( k1=0 ; k1<Dim ; k1++ )
        for ( k2=0 ; k2<Dim ; k2++ ) 
          Stress_bonds[k1][k2] += dr[k1] * dr[k2] * mf ;

    }//for ( j=0 ; j<n_bonds ;
  }// for ( i=0 ; i<ns_loc

}



