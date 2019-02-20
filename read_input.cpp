#include "globals.h"
void read_anneal( void ) ;

void read_input( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[80] ;

  inp = fopen( "ternary.input" , "r" ) ;


  fscanf( inp , "%d %d" , &Nda , &Ndb ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d" , &phiHA , &Nha ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d" , &phiHB , &Nhb ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %d" , &phiHC , &Nhc ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%lf %lf" , &a_squared , &bB_squared ) ; 
  a_squared *= a_squared ;
  bB_squared *= bB_squared ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%d %lf %lf %lf" , &semiflex, &semiflex_req, 
      &semiflex_k, &semiflex_lam ) ;
  fgets( tt , 80 , inp ) ;
 
  CG_ratio = 1.0 ;

  
  // Blank line //
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &phiP ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Rp ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &Xi ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &B_partics ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &eps ) ;
  fgets( tt , 80 , inp ) ;
  
  // Blank line //
  fgets( tt , 80 , inp ) ;


  ///////////////////////////
  // Simulation parameters //
  ///////////////////////////
  fscanf( inp , "%lf" , &rho0 ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf %lf %lf" , &chiAB, &chiAC, &chiBC ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &kappa ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &mu ) ;
  fgets( tt , 80 , inp ) ;

  Diff = ( double* ) calloc( 4 , sizeof( double ) ) ;
  fscanf( inp , "%lf %lf %lf" , &Diff[0] , &Diff[1], &Diff[3] ) ;
  fgets( tt , 80 , inp ) ;

  // Blank line //
  fgets( tt , 80 , inp ) ;


  ///////////////////////
  // Charge parameters //
  ///////////////////////
  fscanf( inp, "%lf", &Lb) ; fgets (tt , 80, inp) ;
  fscanf( inp, "%lf", &qsmear) ; fgets (tt , 80, inp) ;
  fscanf( inp, "%lf", &Zb1) ; fgets (tt , 80, inp) ;
  fscanf( inp, "%lf", &Zb) ; fgets (tt , 80, inp) ;
  fscanf( inp, "%lf", &Zc) ; fgets (tt , 80, inp) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp , "%lf" , &L[i] ) ; 
  fgets( tt , 80 , inp ) ;

  for ( i=0 ; i<Dim ; i++ ) 
    fscanf( inp , "%d" , &Nx[i] ) ; 
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%lf", &delt ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &pmeorder ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &send_buff ) ;
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%d" , &nsteps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &print_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d %d" , &sample_wait , &sample_freq ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &stress_freq ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%d" , &frame_freq ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%d" , &traj_freq ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%d" , &init_flag ) ;
  if ( init_flag >= 1 )
    fscanf(inp, "%d", &interface_bind_steps ) ;

  fgets( tt , 80 , inp ) ;
  
  if ( myrank == 0 ) printf("stress_freq: %d frame: %d traj: %d init_flag: %d\n" , stress_freq, frame_freq, traj_freq, init_flag ) ;
  fclose( inp ) ;

  read_anneal() ;

}
