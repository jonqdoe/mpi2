#include "globals.h"


void bonds() {

  int i,j, k, id1, id2, present_flag ;
  double mdr2, mdr, dr[Dim], delr, mf ;
  double ub_loc = 0.0 , utmp ;


  for ( i=0 ; i<ns_loc ; i++ ) {
    id1 = my_inds[i] ;

    for ( j=0 ; j<n_bonds[id1] ; j++ ) {
      id2 = bonded_to[ id1 ][ j ] ;

      if ( id2 < id1 )
        continue ;

      mdr2 = pbc_mdr2( x[id1], x[id2], dr ) ;

      mdr = sqrt(mdr2) ;

      delr = mdr - bond_eq[ id1 ][j] ;

      ub_loc += delr * delr * bond_coeff[id1][j] / 2.0 ;
      
      mf = bond_coeff[id1][j] * delr / mdr ;
      
      for ( k=0 ; k<Dim ; k++ ) {
        f[id1][k] -= mf * dr[k] ;
        f[id2][k] += mf * dr[k] ;
      }

    }//for ( j=0 ; j<n_bonds ;
  }// for ( i=0 ; i<ns_loc


  MPI_Allreduce( &ub_loc, &utmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;
  Ubond = utmp ;

}




void bond_stress() {

  int i,j, k1, k2, id1, id2 ;
  double mdr2, mdr, dr[Dim], delr, mf ;

  for ( i=0 ; i<ns_loc ; i++ ) {
    id1 = my_inds[i] ;

    for ( j=0 ; j<n_bonds[id1] ; j++ ) {
      id2 = bonded_to[ id1 ][ j ] ;

      mdr2 = pbc_mdr2( x[id1], x[id2], dr ) ;

      mdr = sqrt(mdr2) ;

      delr = mdr - bond_eq[ id1 ][j] ;
      
      mf = bond_coeff[id1][j] * delr / mdr ;
      for ( k1=0 ; k1<Dim ; k1++ )
        for ( k2=0 ; k2<Dim ; k2++ ) 
          Stress_bonds[k1][k2] += dr[k1] * dr[k2] * mf ;

    }//for ( j=0 ; j<n_bonds ;
  }// for ( i=0 ; i<ns_loc

}






/* 
 * bonds_init sets up the bond list. 
 *
 * All particles only keep track of particles they are bonded to
 * with a higher index.
 */

void bonds_init( ) {

  int i, j, k, m, ind, id2 ;

  for ( i=0 ; i<nstot ; i++ ) 
    n_bonds[i] = 0 ;


  // Diblock bonds //
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=0 ; m<Nda + Ndb ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      if ( m < Nda + Ndb - 1 ) {
        id2 = ind + 1 ;
 
        bonded_to[ ind ][ n_bonds[ind] ] = id2 ;
        bond_eq[ ind ][ n_bonds[ind] ] = 0.0 ;
        if ( ind < Nda )
          bond_coeff[ ind ][ n_bonds[ind] ] = 3.0 ;
        else
          bond_coeff[ ind ][ n_bonds[ind] ] = 3.0 ;
        n_bonds[ind]++ ;
 
        bonded_to[ id2 ][ n_bonds[id2] ] = ind ;
        bond_eq[ id2 ][ n_bonds[id2] ] = 0.0 ;
        if ( id2 < Nda )
          bond_coeff[ id2 ][ n_bonds[id2] ] = 3.0 ;
        else
          bond_coeff[ id2 ][ n_bonds[id2] ] = 3.0 ;
        n_bonds[id2]++ ;
      }
      
    }
  } // for ( i=0 ; i<nT[k]


  // Homopolymer A bonds //
  for ( i=0 ; i<nA ; i++ ) {
    for ( m=0 ; m<Nha ; m++ ) {

      ind = nD * (Nda + Ndb) + i * Nha + m ;
      
      if ( m < Nha - 1 ) {
        id2 = ind + 1 ;

        bonded_to[ ind ][ n_bonds[ind] ] = id2 ;
        bond_eq[ ind ][ n_bonds[ind] ] = 0.0 ;
        bond_coeff[ ind ][ n_bonds[ind] ] = 3.0 ;
        n_bonds[ind]++ ;
        
        bonded_to[ id2 ][ n_bonds[id2] ] = ind ;
        bond_eq[ id2 ][ n_bonds[id2] ] = 0.0 ;
        bond_coeff[ id2 ][ n_bonds[id2] ] = 3.0 ;
        n_bonds[id2]++ ;
      }
    }
  } // for ( i=0 ; i<nT[k]


  // Homopolymer B bonds //
  for ( i=0 ; i<nB ; i++ ) {
    for ( m=0 ; m<Nhb ; m++ ) {

      ind = nD * (Nda + Ndb) + nA * Nha + i * Nhb + m ;
      if ( m < Nhb - 1 ) {
        id2 = ind + 1 ;

        bonded_to[ ind ][ n_bonds[ind] ] = id2 ;
        bond_eq[ ind ][ n_bonds[ind] ] = 0.0 ;
        bond_coeff[ ind ][ n_bonds[ind] ] = 3.0 ;
        n_bonds[ind]++ ;
        
        bonded_to[ id2 ][ n_bonds[id2] ] = ind ;
        bond_eq[ id2 ][ n_bonds[id2] ] = 0.0 ;
        bond_coeff[ id2 ][ n_bonds[id2] ] = 3.0 ;
        n_bonds[id2]++ ;
      }
      
    }
  } // for ( i=0 ; i<nT[k]


  // Homopolymer C bonds //
  for ( i=0 ; i<nC ; i++ ) {
    for ( m=0 ; m<Nhc ; m++ ) {

      ind = nD * (Nda + Ndb) + nA * Nha + nB * Nhb + i * Nhc + m ;
      
      if ( m < Nhc - 1 ) {
        id2 = ind + 1 ;

        bonded_to[ ind ][ n_bonds[ind] ] = id2 ;
        bond_eq[ ind ][ n_bonds[ind] ] = 0.0 ;
        bond_coeff[ ind ][ n_bonds[ind] ] = 3.0 ;
        n_bonds[ind]++ ;
        
        bonded_to[ id2 ][ n_bonds[id2] ] = ind ;
        bond_eq[ id2 ][ n_bonds[id2] ] = 0.0 ;
        bond_coeff[ id2 ][ n_bonds[id2] ] = 3.0 ;
        n_bonds[id2]++ ;
      }
    }
  } // for ( i=0 ; i<nT[k]

}




