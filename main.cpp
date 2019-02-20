#define MAIN
#include "globals.h"
#include "pair_style_gaussian.h"

#include "timing.h"

void update_positions( void ) ;
void initialize( void ) ;
void write_lammps_traj( void ) ;
void write_grid( void ) ;
void forces( void ) ;
double integrate( double* ) ;
void write_stress( void ) ;
void bond_stress( void ) ;
void calc_Unb( void ) ;
void anneal_update( void ) ;
void make_eigenval_map( void ) ;




int main( int argc , char** argv ) {

  int i,j,k ;
  
  main_t_in = time(0) ;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank ) ;
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs ) ;
  fftw_mpi_init() ;

  initialize() ;

  write_lammps_traj( ) ;

  FILE *otp ;
  if ( myrank == 0 ) 
    otp = fopen( "data.dat" , "w" ) ;


  if ( myrank == 0 ) {
    printf("Entering main loop!\n") ; fflush( stdout ) ;
  }



  for ( step = 0 ; step < nsteps ; step++ )  {


    if ( do_anneal && step == next_anneal_update )
      anneal_update() ;


    ////////////////////
    // Core algorithm //
    ////////////////////
    forces() ;


    update_positions() ;


    ///////////////////////
    // Write stress data //
    ///////////////////////
    if ( stress_freq > 0 && step % stress_freq == 0 ) {
      bond_stress() ;

      for ( j=0 ; j<Dim ; j++ ) 
        for ( k=0 ; k<Dim ; k++ ) 
          sts_buf[buff_ind][j][k] = Stress_bonds[j][k] ;

      buff_ind++ ;
    }



    ////////////////////////////////
    // Calculate structure factor //
    ////////////////////////////////
    if ( step > sample_wait && step % sample_freq == 0 ) {
      fftw_fwd( rho[0] , ktmp ) ;
      for ( i=0 ; i<ML ; i++ ) {
        double k2, kv[Dim] ;
        k2 = get_k( i , kv ) ;
        avg_sk[0][i] += ktmp[i] * conj(ktmp[i]) ;
      }
      num_averages += 1.0 ;
    }



    /////////////////////////
    // Write lammps output //
    /////////////////////////
    if ( step % traj_freq == 0 && traj_freq > 0 )
      write_lammps_traj() ;






    ///////////////////
    // Write outputs //
    ///////////////////
    if ( step % print_freq == 0 || step == nsteps-1 ) {
      calc_Unb() ;

      Gaussian gauss_AB( chiAB/rho0, 2.0*a_squared, ML, rho[0], rho[1] ) ;
      cout << "Energied! " << gauss_AB.calc_energy() << " " << U_chiab_gg << endl;
      exit(0);

      if ( mu != 0.0 ) 
        make_eigenval_map() ;

      if ( myrank == 0 ) {
        printf("step %d of %d  Ubond: %lf " , step , nsteps , Ubond ) ;
        if ( semiflex )
          printf("Uangle: %lf ", Uangle ) ;
        
        if ( chiAB > 0.0 ) 
          printf("Uchiab: %lf ", U_chiab_gg + U_chi_pg + U_chi_pp) ;
        
        if ( chiBC > 0.0 ) 
          printf("Uchibc: %lf ", U_chibc);
        
        if ( chiAC > 0.0 ) 
          printf("Uchiac: %lf ", U_chiac);

        if ( kappa > 0.0 ) 
          printf("Ukappa: %lf ", U_kappa_gg + U_kappa_pg + U_kappa_pp ) ;
        
        if ( mu != 0.0 )
          printf("U_ms: %lf lc: %lf", U_ms, lc_order_param ) ;
        printf("\n");
        fflush( stdout ) ;
      }

      //cout << "Proc " << myrank << " now has " << ns_loc << endl; 

      if ( stress_freq > 0 )
        write_stress() ;


      write_grid_data( "rhoda" , rhoda ) ;
      write_grid_data( "rhodb" , rhodb ) ;
      
      fftw_fwd( rho[0], ktmp ) ;
      for ( i=0 ; i<ML ; i++ )
        ktmp[i] = ktmp[i] * conj(ktmp[i]) ;
      write_kspace_data( "sk", ktmp ) ;

      if ( mu != 0.0 ) 
        write_Sfield_data( "Sfield", S_field ) ;

      if ( nA > 0.0 )
        write_grid_data( "rhoha" , rhoha ) ;

      if ( nB > 0.0 )
        write_grid_data( "rhohb" , rhohb ) ;

      if ( nC > 0.0 )
        write_grid_data( "rhohc", rhohc ) ;

      if ( nP > 0.0 ) 
        write_grid_data( "rhop" , rhop ) ;

      if ( Lb > 0.0 )
        write_grid_data( "rhoq", rhoq ) ;

      if ( step > sample_wait ) {
        for ( i=0 ; i<ML ; i++ )
          ktmp2[i] = avg_sk[0][i] / num_averages ;
 
        write_kspace_data( "avg_sk" , ktmp2 ) ;
      }


      if ( myrank == 0 ) {
        fprintf( otp , "%d %lf %lf %lf " , step , Ubond , U_chiab_gg, U_kappa_gg ) ;
        if ( chiAC > 0.0 ) 
          fprintf( otp , "%lf ", U_chiac ) ;
        if ( chiBC > 0.0 ) 
          fprintf( otp , "%lf ", U_chibc ) ;
        if ( semiflex )
          fprintf( otp , "%lf ", Uangle) ;
        if ( mu != 0.0 ) 
          fprintf( otp, "%lf %lf ", U_ms, lc_order_param) ;

        fprintf(otp,"\n");
        fflush( otp ) ;
      }
  
    }// if step % print_Freq == 0
    


    ///////////////////////////////////////
    // Logarithmic spaced configurations //
    ///////////////////////////////////////
    if ( frame_freq > 0 && step == frame_freq ) {
      char nm[20] ;
      if ( nA > 0.0 ) {
        sprintf( nm , "rhoha.frame%d" , step ) ;
        write_grid_data( nm , rhoha ) ;
      }
      if ( nB > 0.0 ) {
        sprintf( nm , "rhohb.frame%d" , step ) ;
        write_grid_data( nm , rhohb ) ;
      }
      if ( nP > 0.0 ) {
        sprintf( nm , "rhop.frame%d" , step ) ;
        write_grid_data( nm , rhop ) ;
      }
      if ( nD > 0.0 ) {
        sprintf( nm , "rhoda.frame%d" , step ) ;
        write_grid_data( nm , rhoda ) ;
        sprintf( nm , "rhodb.frame%d" , step ) ;
        write_grid_data( nm , rhodb ) ;
      }

      frame_freq *= 2 ;
    }



  }

  if ( myrank == 0 ) 
    fclose( otp ) ;

  main_t_out = time(0); 
  if ( myrank == 0 ) {
    printf("Total time: %d mins, tot seconds: %ds\n", (main_t_out - main_t_in)/60, 
        (main_t_out-main_t_in) ) ;
 
    printf("FFT time: %d mins %d secs\n", fft_tot_time/60, fft_tot_time%60 ) ;
    printf("Update time: %dmins %d secs\n", move_tot_time/60, move_tot_time%60 ) ;
    printf("Grid time: %dmins %d secs, tot seconds: %ds\n", grid_tot_time/60, grid_tot_time%60, grid_tot_time ) ;

    if ( nprocs > 1 ) {
      printf("\nUpdate time, not including comm: %d mins %d secs\n", move_minus_comm/60, move_minus_comm%60);
      printf("Total time spent in swap/comm routines: %d mins %d secs, tot sec: %d\n", swap_tot_time/60 , swap_tot_time%60, swap_tot_time ) ;
      printf("Total time for exchanging forces: %d s, ghosts: %d s, partics: %d s\n",
          force_comm_tot_time, swap_ghosts_tot_time, swap_partics_tot_time );
      if ( time_debug_tot_time > 0 ) 
        printf("\nDebug total time: %d\n", time_debug_tot_time ) ;
    }
  }

  fftw_mpi_cleanup();

  MPI_Finalize() ;

  return 0 ;

}
