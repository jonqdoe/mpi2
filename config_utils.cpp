#include "globals.h"
void charge_grid( void ) ;
void write_grid_data( const char* , double* ) ;

void random_config( void ) ;
void interface_config( void ) ;
void write_gro( void ) ;
void read_input_conf( FILE* ) ;

void send_all_x_to_root( void ) ;
int mpi_Bcast_posits( int, int, double**, int, MPI_Comm) ;


void initialize_configuration( ) {

  int i, j ;

  if ( myrank == 0 ) {
    if ( init_flag == 0 ) 
      random_config() ;
    else if ( init_flag >= 1 )
      interface_config() ;

    FILE *inp ;
    inp = fopen( "input.lammpstrj" , "r" ) ;
    if ( inp != NULL ) {
      read_input_conf( inp ) ;
      fclose( inp ) ;
      printf("input.lammpstrj read!\n" ) ;
    }
  }

  MPI_Barrier( MPI_COMM_WORLD ) ;

  if ( nprocs > 1 ) {
    mpi_Bcast_posits( nstot, Dim, x, 0, MPI_COMM_WORLD ) ;
    MPI_Bcast( &(tp[0]), nstot, MPI_INT, 0, MPI_COMM_WORLD ) ;
  }
  
  MPI_Barrier( MPI_COMM_WORLD ) ;

}







void read_input_conf( FILE *ip ) {


  int i, j, k, di ;
  double df ;
  int rtflag ;
  char tt[80] ;
  fgets( tt , 80 , ip ) ;
  fgets( tt , 80 , ip ) ;
  fgets( tt , 80 , ip ) ;
  fscanf( ip , "%d\n" , &di ) ;

  if ( di != nstot ) 
    die("Number of sites in input.lammpstrj does not match!\n");

  printf("\nUsing positions from input.lammpstrj!\n\n");
  
  fgets( tt , 80 , ip ) ;

  fgets( tt , 80 , ip ) ;
  fgets( tt , 80 , ip ) ;
  fgets( tt , 80 , ip ) ;
  
  fgets( tt , 80 , ip ) ;

  for ( i=0 ; i<nstot ; i++ ) {
    fscanf( ip, "%d", &di ) ;
    fscanf( ip, "%d", &di ) ;
    fscanf( ip, "%d", &di ) ;
    for ( j=0 ; j<Dim ; j++ ) 
      fscanf( ip , "%lf", &x[i][j] ) ;

    // Read the z position if a 2D simulation
    for ( j=Dim ; j<3 ; j++ )
      fscanf( ip , "%lf" , &df ) ;

    fgets( tt, 80 , ip ) ;
  }

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
    for ( k=0 ; k<nD ; k++ ) {
      for ( m=0 ; m<Nda + Ndb ; m++ ) {
        fprintf( otp , "%d %d %d", ind+1, tp[ind], resind+1 ) ;
   
        for ( j=0 ; j<Dim ; j++ )
          fprintf( otp , " %lf" , x[ind][j]  ) ;
        
        for ( j=Dim ; j<3 ; j++ )
          fprintf( otp , " %lf" , 0.0 );
        
        fprintf( otp , "\n" ) ;
   
        ind++;
      }
      resind++ ;
    }
  
    for ( k=0 ; k<nA ; k++ ) {
      for ( m=0 ; m<Nha ; m++ ) {
        fprintf( otp , "%d %d %d", ind+1, tp[ind], resind+1 ) ;
   
        for ( j=0 ; j<Dim ; j++ )
          fprintf( otp , " %lf" , x[ind][j]  ) ;
        
        for ( j=Dim ; j<3 ; j++ )
          fprintf( otp , " %lf" , 0.0 );
        
        fprintf( otp , "\n" ) ;
   
        ind++;
      }
      resind++ ;
    }
  
    for ( k=0 ; k<nB ; k++ ) {
      for ( m=0 ; m<Nhb ; m++ ) {
        fprintf( otp , "%d %d %d", ind+1, tp[ind], resind+1 ) ;
   
        for ( j=0 ; j<Dim ; j++ )
          fprintf( otp , " %lf" , x[ind][j]  ) ;
        
        for ( j=Dim ; j<3 ; j++ )
          fprintf( otp , " %lf" , 0.0 );
        
        fprintf( otp , "\n" ) ;
   
        ind++;
      }
      resind++ ;
    }
  
    for ( k=0 ; k<nC ; k++ ) {
      for ( m=0 ; m<Nhc ; m++ ) {
        fprintf( otp , "%d %d %d", ind+1, tp[ind], resind+1 ) ;
   
        for ( j=0 ; j<Dim ; j++ )
          fprintf( otp , " %lf" , x[ind][j]  ) ;
        
        for ( j=Dim ; j<3 ; j++ )
          fprintf( otp , " %lf" , 0.0 );
        
        fprintf( otp , "\n" ) ;
   
        ind++;
      }
      resind++ ;
    }
  
    for ( k=0 ; k<nP; k++ ) {
      fprintf( otp , "%d %d %d", ind+1, tp[ind], resind+1 ) ;
   
      for ( j=0 ; j<Dim ; j++ )
        fprintf( otp , " %lf" , x[ind][j]  ) ;
      
      for ( j=Dim ; j<3 ; j++ )
        fprintf( otp , " %lf" , 0.0 );
      
      fprintf( otp , "\n" ) ;
      ind++ ;
      resind++ ;
    }
   
  
    fclose(otp) ;

  }// if ( myrank == 0 ) 

  MPI_Barrier(MPI_COMM_WORLD ) ;
}



void interface_config() {

  int i, m, j, k, ind = 0, n ;

  cout << "\nSETTING UP INTERFACE CONFIGURATION!!" << endl;
  double A = 1.0 ;
  for ( i=1 ; i<Dim; i++)
    A *= L[i] ;
  cout << "Diblock areal density: " << double(nD)/A << " ch/b^2 " 
    << double(nD*(Nda+Ndb-1))/A/6.0 << " ch/Rg^2" << endl << endl;

  // Place the diblocks
  for ( k=0 ; k<nD ; k++ ) { 

    for ( j=0 ; j<Dim ; j++ )
      x[ind][j] = ran2() * L[j] ;
    
    x[ind][0] = 2.0 * ( ran2() - 0.5 ) + L[0]/2.0 ;
    
    // Make two stripes //
    if ( init_flag == 2 ) {
      if ( ran2() < 0.5 )
        x[ind][1] = 0.25 * L[1] * ran2() ;
      else
        x[ind][1] = 0.25 * L[1] * ran2() + 0.5 * L[1] ;
    }

    // Make one stripe //
    if ( init_flag == 5 ) {
      if ( Dim == 3 )
        x[ind][2] = 0.5 * L[2] * ran2() ;
      else
        x[ind][1] = 0.5 * L[1] * ran2() ;
    }
    
    // Make three stripes //
    if ( init_flag == 3 ) {
      double r2 = ran2() ;
      if ( r2 < 0.333 )
        x[ind][1] = 0.167 * L[1] * ran2() ;
      else if ( r2 < 0.666 )
        x[ind][1] = 0.167 * L[1] * ran2() + 0.33 * L[1] ;
      else 
        x[ind][1] = 0.167 * L[1] * ran2() + 0.66 * L[1] ;
    }

    // Make four droplets //
    if ( init_flag == 4 ) {
      if ( Dim != 3 )
        die("Cannot use init_flag == 4 in less than 3D!");

      double r2 = ran2() ;
      if ( r2 < 0.25 ) {
        x[ind][1] = 0.25 * L[1] * ran2() ;
        x[ind][2] = 0.25 * L[2] * ran2() ;
      }
      else if ( r2 < 0.5 ) {
        x[ind][1] = 0.25 * L[1] * ran2() + 0.5 * L[1] ;
        x[ind][2] = 0.25 * L[2] * ran2() ;
      }
      else if ( r2 < 0.75 ) {
        x[ind][1] = 0.25 * L[1] * ran2() ;
        x[ind][2] = 0.25 * L[2] * ran2() + 0.5 * L[2] ;
      }
      else {
        x[ind][1] = 0.25 * L[1] * ran2() + 0.5 * L[1] ;
        x[ind][2] = 0.25 * L[2] * ran2() + 0.5 * L[2] ;
      }
    }

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
    mdu = 1.0 * sqrt( mdu ) ;
    if ( rand_u[0] < 0.0 )
      rand_u[0] *= -1.0 ;

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


  // Random A homopolymers in x < L[0]/2.
  for ( k=0 ; k<nA ; k++ ) {
    for ( j=1 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;
    x[ind][0] = ran2()  * L[0]/2.0 ;

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


  // Random C homopolymers in x > L[0] / 2.
  for ( k=0 ; k<nC ; k++ ) {
    for ( j=0 ; j<Dim ; j++ ) 
      x[ind][j] = ran2() * L[j] ;
    x[ind][0] = ran2() * L[0]/2.0 + L[0]/2.0 ;
    
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
    double min_mdr2 = 0.0, dr[Dim], mdr2 ;

    for ( j=0 ; j<Dim ; j++ )
      x[ind][j] = ran2() * L[j] ;
    x[ind][0] = 2.0 * ( ran2() - 0.5 ) + L[0]/2.0 ;
    
    // Make two stripes //
    if ( init_flag == 2 ) {
      if ( ran2() < 0.5 )
        x[ind][1] = 0.25 * L[1] * ran2() + 0.25 * L[1] ;
      else
        x[ind][1] = 0.25 * L[1] * ran2() + 0.75 * L[1] ;
    }
    
    // Make one stripe //
    if ( init_flag == 5 ) {
      if ( Dim == 3 )
        x[ind][2] = 0.5 * L[2] * ran2() + 0.5 * L[2] ;
      else
        x[ind][1] = 0.5 * L[1] * ran2() + 0.5 * L[1] ;
    }

    tp[ ind ] = 2 ;

    min_mdr2 = 34382.0 ;
    for ( int id2 = ind-k ; id2 < ind ; id2++ ) {
      mdr2 = pbc_mdr2( x[ind], x[id2], dr ) ;
      if ( mdr2 < min_mdr2 )
        min_mdr2 = mdr2 ;
    }

    if ( k == 0 )
      ind += 1 ;
    else if ( min_mdr2 > 2.0 * Rp )
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
}// interface_config() 


void random_config( void ) {

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
}

