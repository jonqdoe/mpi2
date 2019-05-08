#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std ;

#define PI   3.141592653589793238462643383




int main(int argc, char** argv ) {

  int i, m, j, k, ind = 0, n ;


  for ( k=0 ; k<nD ; k++ ) { 
    for ( j=0 ; j<Dim ; j++ )
      x[ind][j] = ran2() * L[j] ;
    tp[ind] = ( Nda > 0 ? 0 : 1 ) ;

    ind++ ;
      
    for ( m=1 ; m<Nda ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;
 
        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }
 
      tp[ ind ] = ( m < Nda ? 0 : 1 ) ;
 
      ind++ ;
      
    }
    double rand_u[Dim], mdu = 0.0 ;
    for ( j=0 ; j<Dim ; j++ ) {
      rand_u[j] = gasdev2() ;
      mdu += rand_u[j] * rand_u[j] ;
    }
    mdu = sqrt( mdu ) ;

    for ( m=Nda ; m<Nda + Ndb ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + semiflex_req * rand_u[j] / mdu ;
 
        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }
 
      tp[ ind ] = ( m < Nda ? 0 : 1 ) ;
 
      ind++ ;
      
    }
  }// for n=0:Nd-1


  // Random A homopolymers //
  for ( k=0 ; k<nA ; k++ ) {
    for ( j=0 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;

    tp[ ind ] = 0 ;

    ind += 1 ;
 
    for ( m=1 ; m<Nha ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }

      tp[ ind ] = 0 ;

      ind++ ;
    
    }
  }


  // Random B homopolymers //
  for ( k=0 ; k<nB ; k++ ) {
    for ( j=0 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;

    tp[ ind ] = 1 ;

    ind += 1 ;
 
    for ( m=1 ; m<Nhb ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }

      tp[ ind ] = 1 ;

      ind++ ;
    
    }
  }


  // Random C homopolymers //
  for ( k=0 ; k<nC ; k++ ) {
    for ( j=0 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;
    
    tp[ ind ] = 3 ;

    ind += 1 ;
 
    for ( m=1 ; m<Nhc ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;

        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }

      tp[ ind ] = 3 ;

      ind++ ;
    
    }
  }

  // Random particle centers //
  for ( k=0 ; k<nP ; k++ ) {
    for ( j=0 ; j<Dim ; j++ )
      x[ind][j] = ran2() * L[j] ;

    tp[ ind ] = 2 ;

    ind += 1 ;
  }


  // Assign the labels //
  for ( i=0 ; i<nstot ; i++ ) {
    if ( tp[i] == 0 ) 
      xc[i] = "H" ;
    else if ( tp[i] == 1 )
      xc[i] = "He" ;
    else if ( tp[i] == 2 )
      xc[i] = "O" ;
    else if ( tp[i] == 3 )
      xc[i] = "S" ;
    else if ( tp[i] == 4 )
      xc[i] = "N" ;
    else if ( tp[i] == 5 )
      xc[i] = "Br" ;
    else if ( tp[i] == 6 )
      xc[i] = "C" ;
    else if ( tp[i] == 7 )
      xc[i] = "Na" ;
    else if ( tp[i] == 8 )
      xc[i] = "P" ;
    else if ( tp[i] == 9 )
      xc[i] = "Ca" ;
  }
  printf("Random config generated!\n") ; fflush( stdout ) ;
  return 0 ;

}

