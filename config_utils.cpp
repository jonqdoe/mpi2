#include "globals.h"
void allocate_particles( void ) ;
void send_all_x_to_root( void ) ;

void read_config() {

  int i, di, ind, ltp;
  double dx, dy, dz;
  char tt[120];

  FILE *inp;
  inp = fopen("input.data" , "r" );
  if ( inp == NULL ) die("Failed to find input.data!" );
  
  fgets( tt, 120 , inp );
  fgets( tt, 120 , inp );

  fscanf( inp , "%d" , &nstot );      fgets( tt , 120 , inp );
  fscanf( inp , "%d" , &n_total_bonds );  fgets( tt , 120 , inp );
  fscanf( inp , "%d" , &n_total_angles );  fgets( tt , 120 , inp );

  fgets( tt , 120 , inp );

  fscanf( inp , "%d" , &ntypes );  fgets( tt , 120 , inp );
  fscanf( inp , "%d" , &nbond_types );  fgets( tt , 120 , inp );
  fscanf( inp , "%d" , &nangle_types );  fgets( tt , 120 , inp );

  // Read in box shape
  fgets(tt, 120 , inp );
  fscanf( inp , "%lf %lf" , &dx , &dy ) ;   fgets( tt, 120 , inp );
  L[0] = (double) ( dy - dx ) ;
  fscanf( inp , "%lf %lf" , &dx , &dy ) ;   fgets( tt, 120 , inp );
  L[1] = (double) ( dy - dx ) ;
  fscanf( inp , "%lf %lf" , &dx , &dy ) ;   fgets( tt, 120 , inp );
  if ( Dim > 2 )
    L[2] = (double) ( dy - dx ) ;

  V = 1.0 ;
  for ( i=0 ; i<Dim ; i++ ) {
    Lh[i] = 0.5 * L[i] ;
    V *= L[i] ;
  }
  
  
  // Allocate memory for positions //
  allocate_particles();
  if ( myrank == 0 )
    printf("Particle memory allocated!\n");

  fgets(tt, 120 , inp );

  // Read in particle masses
  fgets(tt, 120 , inp );
  fgets(tt, 120 , inp );
  for ( i=0 ; i<ntypes ; i++ ) {
    fscanf( inp , "%d %lf" , &di , &dx ); fgets( tt, 120, inp );
    mass[ di-1 ] = ( double ) dx ;
  }
  fgets( tt, 120, inp );


  // Read in atomic positions
  fgets( tt, 120 , inp );
  fgets( tt, 120 , inp );


  for ( i=0 ; i<nstot ; i++ ) {
    if ( feof(inp) ) die("Premature end of input.conf!");

    fscanf( inp , "%d %d %d", &ind , &di , &ltp );
    ind -= 1;

    molecID[ind] = di-1 ;
    tp[ind] = ltp-1;

    for (int  j=0 ; j<Dim ; j++ ) {
      fscanf( inp, "%lf", &dx ) ;
      x[ind][j] = dx ;
    }

    fgets( tt , 120 , inp );
    if ( i == nstot-1 )
      printf("%lf %lf\n", x[ind][0], x[ind][1] ) ;
  }
  fgets( tt , 120 , inp );


  // Read in bond information
  fgets( tt , 120 , inp );
  fgets( tt , 120 , inp );

  for ( i=0 ; i<nstot ; i++ )
    n_bonds[i] = 0 ;

  for ( i=0 ; i<n_total_bonds ; i++ ) {
    fscanf( inp , "%d" , &di );
    fscanf( inp , "%d" , &di );
    int b_type = di-1;
    
    fscanf( inp , "%d" , &di );
    int i1 = di-1;
    
    fscanf( inp , "%d" , &di );
    int i2 = di - 1;

    if ( i2 < i1 ) {
      di = i2 ;
      i2 = i1 ;
      i1 = di ;
    }

    bonded_to[i1][ n_bonds[i1] ] = i2 ;
    bond_type[i1][ n_bonds[i1] ] = b_type ;
    n_bonds[i1]++ ;

    bonded_to[i2][ n_bonds[i2] ] = i1 ;
    bond_type[i2][ n_bonds[i2] ] = b_type ;
    n_bonds[i2]++ ;

  }
  fgets( tt , 120 , inp );



  // Read in angle information
  fgets( tt , 120 , inp );
  fgets( tt , 120 , inp );
  for ( i=0 ; i<n_total_angles ; i++ ) {
    
    fscanf( inp , "%d" , &di );
    fscanf( inp , "%d" , &di );
    
    int a_type = di - 1;

    fscanf( inp , "%d" , &di );
    int i1 = di - 1;

    fscanf( inp , "%d" , &di );
    int i2 = di - 1;

    fscanf( inp , "%d" , &di );
    int i3 = di - 1;

    if ( i3 < i1 ) {
      di = i3 ;
      i3 = i1 ;
      i1 = di ;
    }

    int na = n_angles[i1] ;
    angle_first[i1][na] = i1 ;
    angle_mid[i1][na] = i2 ;
    angle_end[i1][na] = i3 ;
    angle_type[i1][na] = a_type ;
    n_angles[i1] += 1 ;

    na = n_angles[i2] ;
    angle_first[i2][na] = i1 ;
    angle_mid[i2][na] = i2 ;
    angle_end[i2][na] = i3 ;
    angle_type[i2][na] = a_type ;
    n_angles[i2] += 1 ;

    na = n_angles[i3] ;
    angle_first[i3][na] = i1 ;
    angle_mid[i3][na] = i2 ;
    angle_end[i3][na] = i3 ;
    angle_type[i3][na] = a_type ;
    n_angles[i3] += 1 ;

    fgets( tt , 120 , inp );

  }
 
  fclose( inp );

}








void write_lammps_traj() {

  send_all_x_to_root() ;

  if ( myrank == 0 ) {
    FILE *otp ;
    int i, j, ind, k, m, n, resind ;
    if ( step == 0 )
      otp = fopen( "traj.lammpstrj", "w" ) ;
    else
      otp = fopen( "traj.lammpstrj", "a" ) ;
  
    fprintf( otp , "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n", step, nstot ) ;
    fprintf( otp , "ITEM: BOX BOUNDS pp pp pp\n" ) ;
    fprintf( otp , "0.0 %lf\n0.0 %lf\n0.0 %lf\n", L[0], L[1], (Dim==3 ? L[2] : 1.0 ) ) ;
  
    fprintf( otp , "ITEM: ATOMS id type mol x y z\n" ) ;
    ind = 0 ;
    resind = 0 ;
    for ( i=0 ; i<nstot; i++ ) {
      fprintf( otp, "%d %d %d  ", i+1, tp[i], molecID[i] ) ;
      for ( j=0 ; j<Dim ; j++ ) 
        fprintf(otp, "%lf ", x[i][j] ) ;

      for ( j=Dim ; j<3 ; j++ ) 
        fprintf(otp, "%lf", x[i][j] ) ;
      fprintf(otp, "\n");

    }
  
    fclose(otp) ;

  }// if ( myrank == 0 ) 

  MPI_Barrier(MPI_COMM_WORLD ) ;
}




