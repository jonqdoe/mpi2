#include "globals.h"
#include "timing.h"

void allocate_grid( void ) ;
int fft_init( void ) ;
void charge_grid( void ) ;
void read_input( void ) ;
void swap_ghosts( void ) ;
void find_my_particles( void ) ;
void read_config( void ) ;
void initialize_potential( void ) ;



void initialize() {

  if ( nprocs == 2 ) 
    die("Invalid number of processors! This stupid code isn't set up to deal with only two procs!");

  int i,j ;

  idum =  -long( time(0) ) * ( myrank + 1 ) ; // 10 ; //

  for ( i=0 ; i<nprocs ; i++ ) 
    if ( i == myrank ) {
      printf("rank: %d idum: %ld\n", myrank, idum) ;
      fflush(stdout);
    }

  move_minus_comm = move_tot_time = fft_tot_time = grid_tot_time = 0 ;
  force_comm_tot_time = swap_tot_time = swap_build_list_tot_time = 0 ;
  swap_partics_tot_time = swap_ghosts_tot_time = time_debug_tot_time = 0 ;
 
  read_config() ;

  read_input() ;


  mem_use = 0. ;
  do_anneal = 0 ;

  M = 1 ;
  ML = 1 ;
  grid_per_partic = 1 ;
  for ( j=0 ; j<Dim ; j++ ) {
    M *= Nx[j] ;
    ML *= NxL[j] ;
    dx[j] = L[j] / double( Nx[j] ) ;

    grid_per_partic *= ( pmeorder + 1 ) ;
    if ( myrank == 0 ) 
      printf("dx: %lf\n" , dx[j] ) ;
  }
  gvol = V / double( M ) ;

  int fft_alloc = fft_init() ;
  if ( myrank == 0 ) printf("FFTW-MPI Initialized\n") ; fflush( stdout ) ;

  z_min = double(zstart) * dx[Dim-1] ;
  z_max = double(zstart + NxL[Dim-1]) * dx[Dim-1] ;
  printf("Z range for process %d: [%lf, %lf)\n", myrank,z_min,z_max ) ;
  MPI_Barrier(MPI_COMM_WORLD);

  if ( z_max - z_min < 2.0 * buff_size )
    die("OVERLAPPING COMMUNICATION REGIONS!\nUse fewer processors or run a more dense grid.\n") ;


  step = 0 ;
  
  buff_ind = 0 ;
  if ( stress_freq > print_freq )
    stress_freq = print_freq ;

  if ( stress_freq > 0 )
    buff_size = print_freq / stress_freq + 1;
  else
    buff_size = 0 ;

  I = complex<double>( 0.0 , 1.0 ) ;

  if ( myrank == 0 ) 
    printf("Total segments: %d\n" , nstot ) ;
  if ( myrank == 0 ) 
    printf("grid_vol: %lf\n" , gvol ) ;
  if ( myrank == 0 ) 
    printf("Particles per grid point: %lf\n" , double(nstot) / double(M) ) ;

  allocate_grid() ;
  
  if ( myrank == 0 ) printf("Memory allocated: %lf MB\n" , mem_use / 1.0E6 ) ; 
  

  find_my_particles() ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) { cout << "Found my particles!\n" ; fflush(stdout) ; }


  printf("0: n_bonds[0]: %d [1]: %d\n", n_bonds[0], n_bonds[1]) ;

  swap_ghosts() ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) { cout << "Exchanged ghost particles!\n" ; fflush(stdout) ; }


  printf("0: n_bonds[0]: %d [1]: %d\n", n_bonds[0], n_bonds[1]) ;

  charge_grid() ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) printf("grid charged\n"); fflush(stdout); 


  printf("0: n_bonds[0]: %d [1]: %d\n", n_bonds[0], n_bonds[1]) ;

  initialize_potential() ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) printf("potentials initialized, written\n") ; fflush( stdout ) ; 

}



void initialize_potential( ) {

  Gauss = new Gaussian[n_gaussian_pairstyles] ;
  for ( int i=0 ; i<n_gaussian_pairstyles ; i++ ) {
    Gauss[i].Initialize_Gaussian( gaussian_prefactor[i], 
        gaussian_sigma[i] * gaussian_sigma[i], 
        ML,
        Components[ gaussian_types[i][0] ], 
        Components[ gaussian_types[i][1] ] ) ;
  }

}



int allocate_mpi_tmp( int, int ) ;

void allocate_particles() {

  int i, j ;
  int rtflag ;

  mem_use += allocate_mpi_tmp( nstot*Dim, ML ) ;

  rtflag = malloc2ddouble( &x, nstot, Dim ) ;
  rtflag = malloc2ddouble( &rec_N_x, nstot, Dim ) ;
  rtflag = malloc2ddouble( &rec_S_x, nstot, Dim ) ;
  rtflag = malloc2ddouble( &send_N_x, nstot, Dim ) ;
  rtflag = malloc2ddouble( &send_S_x, nstot, Dim ) ;

  xc = ( char** ) calloc( nstot , sizeof( char* ) ) ;
  f = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  for ( i=0 ; i<nstot ; i++ ) {
    f[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
    xc[i] = ( char* ) calloc( 8 , sizeof( char ) ) ;
  }

  mem_use += nstot * 2 * 3 * Dim * sizeof( double ) ; 
  mem_use += nstot * 8 * sizeof( char ) ; 
  
  mass =        ( double* ) calloc( nstot, sizeof( double ) ) ;
  molecID =     ( int* ) calloc( nstot, sizeof( int ) ) ;
  tp =          ( int* ) calloc( nstot , sizeof( int ) );
  local_flag =  ( int* ) calloc( nstot , sizeof( int ) );
  my_inds =     ( int* ) calloc( nstot , sizeof( int ) );
  send_N_inds = ( int* ) calloc( nstot , sizeof( int ) );
  send_S_inds = ( int* ) calloc( nstot , sizeof( int ) );
  rec_N_inds =  ( int* ) calloc( nstot , sizeof( int ) );
  rec_S_inds =  ( int* ) calloc( nstot , sizeof( int ) );
  ghost_inds =  ( int* ) calloc( nstot , sizeof( int ) );

  grid_inds = ( int** ) calloc( nstot , sizeof( int* ) );
  grid_W = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  for ( i=0 ; i<nstot ; i++ ) {
    grid_inds[i] = ( int* ) calloc( grid_per_partic , sizeof( int ) ) ;
    grid_W[i] = ( double* ) calloc( grid_per_partic , sizeof( double ) ) ;
  }

  mem_use += ( 2*nstot + nstot * grid_per_partic ) * sizeof( int ) ;
  mem_use += nstot * grid_per_partic * sizeof( double ) ;


  // NOTE: Assumes that a particle is bonded to a maximum of three
  // particles
  n_bonds =     ( int* ) calloc( nstot, sizeof( int ) ) ;
  n_angles =    ( int* ) calloc( nstot, sizeof( int ) ) ;
  bonded_to =   ( int** ) calloc( nstot, sizeof( int* ) ) ;
  bond_type =   ( int** ) calloc( nstot, sizeof( int* ) ) ;
  angle_first = ( int** ) calloc( nstot, sizeof( int* ) ) ;
  angle_mid =   ( int** ) calloc( nstot, sizeof( int* ) ) ;
  angle_end =   ( int** ) calloc( nstot, sizeof( int* ) ) ;
  angle_type =   ( int** ) calloc( nstot, sizeof( int* ) ) ;

  for ( i=0 ; i<nstot ; i++ ) {
    bonded_to[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    bond_type[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    
    angle_first[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    angle_mid[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    angle_end[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    angle_type[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
  }

  mem_use += nstot * ( 5 * sizeof( int ) ) ;


  
  bond_coeff =  ( double* ) calloc( nbond_types, sizeof( double ) ) ;
  bond_eq =     ( double* ) calloc( nbond_types, sizeof( double ) ) ;
  mem_use += nbond_types * 2 * sizeof(double) ;


  angle_coeff = ( double** ) calloc( nangle_types, sizeof( double* ) ) ;
  for ( i=0 ; i<nangle_types ; i++ ) 
    angle_coeff[i] = ( double* ) calloc( 5, sizeof(double) ) ;

  mem_use += nangle_types * 5 * sizeof(double) ;
}


void allocate_grid( ) {

  int i, j ;
  int rtflag ;

  ktmp = ( complex<double>* ) calloc( ML , sizeof( complex<double> ) ) ;
  ktmp2 = ( complex<double>* ) calloc( ML , sizeof( complex<double> ) ) ;
  tmp = ( double* ) calloc( ML , sizeof( double ) ) ;
  tmp2 = ( double* ) calloc( ML , sizeof( double ) ) ;
  
  mem_use += 6 * ML * sizeof( double ) ; 

  Components = new FieldComponent[ntypes] ;
  for ( i=0 ; i<ntypes ; i++ ) 
    Components[i].Initialize( ML ) ;
}

void allocate_gaussians() {
  gaussian_prefactor = ( double* ) calloc( n_gaussian_pairstyles, sizeof(double) ) ;
  gaussian_sigma = ( double* ) calloc( n_gaussian_pairstyles, sizeof(double) ) ;
  gaussian_types = ( int** ) calloc( n_gaussian_pairstyles, sizeof( int* ) ) ;
  for (int i=0 ; i<n_gaussian_pairstyles ; i++ ) 
    gaussian_types[i] = (int*) calloc( 2, sizeof(int) ) ;
}
