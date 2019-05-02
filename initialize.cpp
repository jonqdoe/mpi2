#include "globals.h"
#include "timing.h"

void allocate( void ) ;
int fft_init( void ) ;
void initialize_potential( void ) ;
void initialize_configuration( void ) ;
void charge_grid( void ) ;
void read_input( void ) ;
void swap_ghosts( void ) ;
void find_my_particles( void ) ;
void bonds_init( void ) ;
void angle_init( void ) ;


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

  read_input() ;
  
  if ( phiP + phiHA + phiHB + phiHC > 1.0 )
    die("Invalid volume fractions!\n") ;


  mem_use = 0. ;

  M = 1 ;
  V = 1.0 ;
  ML = 1 ;
  grid_per_partic = 1 ;
  for ( j=0 ; j<Dim ; j++ ) {
    V *= L[j] ;
    M *= Nx[j] ;
    ML *= NxL[j] ;
    dx[j] = L[j] / double( Nx[j] ) ;
    Lh[j] = 0.5 * L[j] ;
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

  // This is used for density fields //
  ntypes = 4 ;

  Rg = sqrt( double( Nda + Ndb ) / 6.0 ) ;
  Rg3 = Rg * Rg * Rg ;
  //rho0 = C * double( Nda + Ndb ) / Rg3 ;


  nD = int( ( 1.0 - phiHA - phiHB - phiHC - phiP ) * rho0 * V / ( Nda + Ndb ) ) ;
  nA = int( phiHA * rho0 * V / Nha ) ;
  nB = int( phiHB * rho0 * V / Nhb ) ;
  nC = int( phiHC * rho0 * V / Nhc ) ;



  Vp = rho0 ;
  if ( Dim == 2 )
    Vp *= PI * Rp * Rp ;
  else if ( Dim == 3 )
    Vp *= 4.0 * PI * Rp * Rp * Rp / 3.0 ;

  Diff[2] = 1.0 / Vp ;
  

  nP = int( phiP * rho0 * V / Vp ) ;

  nsD = nD * (Nda + Ndb) ;
  nsA = nA * Nha ;
  nsB = nB * Nhb ;
  nsC = nC * Nhc ;

  if ( myrank == 0 ) 
    printf("Input rho0: %lf , " , rho0 ) ;

  rho0 = ( nD * (Nda + Ndb ) + nA * Nha + nB * Nhb + nC * Nhc + nP * Vp ) / V  ;
  if ( myrank == 0 ) 
    printf("actual rho0: %lf\n" , rho0 ) ;
  if ( myrank == 0 ) 
    printf("Diblock C: %lf\n", rho0/(Nda+Ndb)*pow(double(Nda+Ndb)/6.0, 1.5) ) ;
  if ( myrank == 0 ) 
    printf("\nnD: %d\nnA: %d\nnB: %d\nnC: %d\nnP: %d\n\n" , nD, nA, nB, nC, nP ) ;

  // Derived quantities //
  nstot = nA * Nha + nB * Nhb + nD * ( Nda + Ndb ) + nC * Nhc + nP ;
 
  step = 0 ;
  num_averages = 0.0 ;
  
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

  allocate() ;
  
  if ( myrank == 0 ) printf("Memory allocated: %lf MB\n" , mem_use / 1.0E6 ) ; 
  
  initialize_configuration() ;
  
  if ( myrank == 0 ) printf("Initial config generated\n") ;


  bonds_init() ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) printf("bonds initialized\n") ;

  if ( semiflex ) {
    angle_init() ;
  
    MPI_Barrier(MPI_COMM_WORLD) ;
    if ( myrank == 0 ) printf("angles initialized\n") ;
  }

  find_my_particles() ;

  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) { cout << "Found my particles!\n" ; fflush(stdout) ; }

  swap_ghosts() ;
  
  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) { cout << "Exchanged ghost particles!\n" ; fflush(stdout) ; }


  charge_grid() ;

  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) printf("grid charged\n"); fflush(stdout); 

  initialize_potential() ;

  MPI_Barrier(MPI_COMM_WORLD) ;
  if ( myrank == 0 ) printf("potentials initialized, written\n") ; fflush( stdout ) ; 

}



void initialize_potential( ) {

  int i, j;

  double ro[Dim], rc[Dim], dr[Dim], mdr2 , pref , mdr , k2, kv[Dim] ;
  pref = V / ( pow( 2.0 * sqrt(PI) , Dim ) ) ; // Note: the factor of V comes from the FFT

  for ( j=0 ; j<Dim ; j++ )
    ro[j] = 0.0 ;

  for ( i=0 ; i<ML ; i++ ) {

    get_r( i , rc ) ;

    mdr2 = pbc_mdr2( ro, rc, dr ) ;
    mdr = sqrt( mdr2 ) ;

    uG[i] = exp( -mdr2 / 4.0 / a_squared ) * pref ;
    r_dudr[i] = -mdr2 * exp( -mdr2 / 2.0 ) ;

    tmp[i] = rho0 / 2.0 * ( 1.0 - erf( ( mdr - Rp ) / Xi ) ) * V;
    gammaP[i] = tmp[i] ;

    if (eps != 0.0 )
      gamma_sig[i] = V * rho0 * exp( -(mdr-Rp-2.0)*(mdr-Rp-2.0)/Xi ) ;
  }


  // Set up the particle-particle potential //
  fftw_fwd( tmp , ktmp ) ;
  for ( i=0 ; i<ML ; i++ ) 
    ktmp2[i] = ktmp[i] * ktmp[i] ;
  fftw_back( ktmp2 , uP ) ;

  // Set up particle-polymer potential //
  for ( i=0 ; i<ML ; i++ ) {
    k2 = get_k( i , kv ) ;
    ktmp[i] *= exp( -k2 / 2.0 ) ;
  }
  fftw_back( ktmp , uPG ) ;

  // Adsorbing polymer/particle interaction
  if ( eps != 0.0 ) {
    fftw_fwd( gamma_sig, ktmp ) ;
    for ( i=0 ; i<ML ; i++ ) {
      k2 = get_k( i, kv ) ;
      ktmp[i] *= exp( -k2 / 2.0 ) ;
    }
    fftw_back( ktmp, uAG ) ;
  }


  for ( j=0 ; j<Dim ; j++ ) {
    field_gradient( uG , grad_uG[j] , j ) ;
    field_gradient( uP , grad_uP[j] , j ) ;
    field_gradient( uPG , grad_uPG[j] , j ) ;
    if ( eps != 0.0 )
      field_gradient( uAG, grad_uAG[j], j ) ;
  }


//  write_grid_data( "ug" , uG ) ;
//  write_grid_data( "up" , uP ) ;
//  write_grid_data( "upg" , uPG ) ;


  for ( j=0 ; j<Dim ; j++ ) {
    char nm[20] ;
    fftw_fwd( grad_uG[j] , grad_uG_hat[j] ) ;
//    sprintf( nm , "grad_ug_%d" , j ) ;
//    write_grid_data( nm , grad_uG[j] ) ;

    fftw_fwd( grad_uP[j] , grad_uP_hat[j] ) ;

    fftw_fwd( grad_uPG[j] , grad_uPG_hat[j] ) ;

    if ( eps != 0.0 ) 
      fftw_fwd( grad_uAG[j], grad_uAG_hat[j] ) ;
  }
}



int allocate_mpi_tmp( int, int ) ;

void allocate( ) {

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

  if ( mu != 0.0 ) {
    mono_u = ( double** ) calloc( nstot, sizeof( double* ) ) ;
    mono_S = ( double*** ) calloc( nstot, sizeof( double** ) ) ;
    for ( i=0 ; i<nstot ; i++ ) {
      mono_u[i] = ( double* ) calloc( Dim, sizeof(double) ) ;
      mono_S[i] = ( double** ) calloc( Dim, sizeof(double*) ) ;
      for ( j=0 ; j<Dim ; j++ ) {
        mono_u[i][j] = 0.0 ;
        mono_S[i][j] = ( double* ) calloc( Dim, sizeof( double ) ) ;
      }
    }
    mem_use += (Dim+1)*Dim * nstot * sizeof(double) ;

    S_field    = ( double*** ) calloc( ML, sizeof( double** ) ) ;
    S_conv_u   = ( double*** ) calloc( ML, sizeof( double** ) ) ;
    tmp_tensor = ( double*** ) calloc( ML, sizeof( double** ) ) ;
    for ( i=0 ; i<ML ; i++ ) {
      S_field[i]    = ( double** ) calloc( Dim, sizeof( double* ) ) ;
      S_conv_u[i]   = ( double** ) calloc( Dim, sizeof( double* ) ) ;
      tmp_tensor[i] = ( double** ) calloc( Dim, sizeof( double* ) ) ;
      for ( j=0 ; j<Dim ; j++ ) {
        S_field[i][j]    = ( double* ) calloc( Dim, sizeof(double) ) ;
        S_conv_u[i][j]   = ( double* ) calloc( Dim, sizeof(double) ) ;
        tmp_tensor[i][j] = ( double* ) calloc( Dim, sizeof(double) ) ;
      }
    }

    mem_use += 3 * Dim * Dim * ML * sizeof( double ) ;

    if ( myrank == 0 )
      printf("Allocated memory for orientation parameters!\n");
  }

  sts_buf = ( double*** ) calloc( buff_size , sizeof( double** ) ) ;
  for ( i=0 ; i<buff_size ; i++ ) {
    sts_buf[i] = ( double** ) calloc( Dim , sizeof( double* ) ) ;
    for ( j=0 ; j<Dim ; j++ )
      sts_buf[i][j] = ( double* ) calloc( Dim , sizeof( double ) ) ;
  }
  mem_use += Dim * Dim * buff_size * sizeof( double ) ;


  // NOTE: Assumes that a particle is bonded to a maximum of three
  // particles
  n_bonds =     ( int* ) calloc( nstot, sizeof( int ) ) ;
  n_angles =    ( int* ) calloc( nstot, sizeof( int ) ) ;
  bonded_to =   ( int** ) calloc( nstot, sizeof( int* ) ) ;
  angle_first = ( int** ) calloc( nstot, sizeof( int* ) ) ;
  angle_mid =   ( int** ) calloc( nstot, sizeof( int* ) ) ;
  angle_end =   ( int** ) calloc( nstot, sizeof( int* ) ) ;
  bond_coeff =  ( double** ) calloc( nstot, sizeof( double* ) ) ;
  angle_coeff = ( double** ) calloc( nstot, sizeof( double* ) ) ;
  bond_eq =     ( double** ) calloc( nstot, sizeof( double* ) ) ;
  for ( i=0 ; i<nstot ; i++ ) {
    bonded_to[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    bond_coeff[i] = ( double* ) calloc( 5, sizeof( double ) ) ;
    bond_eq[i] = ( double* ) calloc( 5, sizeof( double ) ) ;
    
    angle_first[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    angle_mid[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    angle_end[i] = ( int* ) calloc( 5, sizeof( int ) ) ;
    angle_coeff[i] = ( double* ) calloc( 5, sizeof( double ) ) ;
  }

  mem_use += nstot * ( 6 * sizeof( int ) + 10 * sizeof( double ) ) ;



  for ( i=0 ; i<Dim ; i++ ) {
    grad_uG_hat[i] = ( complex<double>* ) calloc( ML , sizeof( complex<double> ) ) ;
    grad_uP_hat[i] = ( complex<double>* ) calloc( ML , sizeof( complex<double> ) ) ;
    grad_uPG_hat[i] = ( complex<double>* ) calloc( ML , sizeof( complex<double> ) ) ;
  }
  
  mem_use += 3 * Dim * ML * sizeof( double ) ; 
  
  if (eps != 0.0 ) {
    uAG = ( double* ) calloc( ML, sizeof( double ) ) ;
    gamma_sig = ( double* ) calloc( ML, sizeof( double ) ) ;
    for ( j=0 ; j<Dim ; j++ ) {
      grad_uAG[j] = ( double* ) calloc( ML, sizeof( double ) ) ;
      grad_uAG_hat[j] = ( complex<double>* ) calloc( ML, sizeof( complex<double> ) ) ;
    }
    mem_use += (Dim + 2 ) * ML * sizeof( double ) ;
    mem_use += Dim * ML * sizeof( complex<double> ) ;
  }


  ktmp = ( complex<double>* ) calloc( ML , sizeof( complex<double> ) ) ;
  ktmp2 = ( complex<double>* ) calloc( ML , sizeof( complex<double> ) ) ;
  tmp = ( double* ) calloc( ML , sizeof( double ) ) ;
  tmp2 = ( double* ) calloc( ML , sizeof( double ) ) ;
  uG = ( double* ) calloc( ML , sizeof( double ) ) ;
  uP = ( double* ) calloc( ML , sizeof( double ) ) ;
  uPG = ( double* ) calloc( ML , sizeof( double ) ) ;
  r_dudr = ( double* ) calloc( ML , sizeof( double ) ) ;
  
  mem_use += 7 * ML * sizeof( double ) ; 

  phiA = ( double* ) calloc( ML , sizeof( double ) ) ;
  Dweights = ( double* ) calloc( nstot, sizeof( double ) ) ;
  mem_use += ML * sizeof( double ) ; 

  rhot = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhoha = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhohb = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhohc = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhoda = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhodb = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhop = ( double* ) calloc( ML , sizeof( double ) ) ;
  gammaP = ( double* ) calloc( ML , sizeof( double ) ) ;
  smrhop = ( double* ) calloc( ML , sizeof( double ) ) ;

  mem_use += 8 * ML * sizeof( double ) ; 
  


  

  rho = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  avg_sk = ( complex<double>** ) calloc( ntypes , sizeof( complex<double>* ) ) ;
  w = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  
  for ( i=0 ; i<ntypes ; i++ ) {
    rho[i] = ( double * ) calloc( ML , sizeof( double ) ) ;
    avg_sk[i] = ( complex<double> * ) calloc( ML , sizeof( complex<double> ) ) ;
    w[i] = ( double* ) calloc( ML , sizeof( double ) ) ;
  }

  mem_use += 2 * ntypes * ML * sizeof( double ) ;
  mem_use += ntypes * ML * sizeof( complex<double> ) ;

  for ( j=0 ; j<Dim ; j++ ) {
    grad_uG[j] = ( double* ) calloc( ML , sizeof( double ) ) ;
    grad_uP[j] = ( double* ) calloc( ML , sizeof( double ) ) ;
    grad_uPG[j] = ( double* ) calloc( ML , sizeof( double ) ) ;
    
    gradwA[j] = ( double* ) calloc( ML , sizeof( double ) ) ;
    gradwB[j] = ( double* ) calloc( ML , sizeof( double ) ) ;
    gradwC[j] = ( double* ) calloc( ML , sizeof( double ) ) ;
    gradwP[j] = ( double* ) calloc( ML , sizeof( double ) ) ;
    
  }

  mem_use += Dim * ML * 8 * sizeof( double ) ;

}
