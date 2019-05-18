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
    cout <<"Entering forces\n" ;
    forces() ;


    update_positions() ;



    ////////////////////////////////
    // Calculate structure factor //
    ////////////////////////////////
    if ( step > sample_wait && step % sample_freq == 0 ) {
      //fftw_fwd( rho[0] , ktmp ) ;
      for ( i=0 ; i<ML ; i++ ) {
        double k2, kv[Dim] ;
        k2 = get_k( i , kv ) ;
        //avg_sk[0][i] += ktmp[i] * conj(ktmp[i]) ;
      }
      //num_averages += 1.0 ;
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
      forces() ;
      calc_Unb() ;

//      printf("making field components...\n"); fflush(stdout) ;
//      FieldComponent A(ML), B(ML) ;
//      FieldComponent Fields[5] ;
//      printf("done!\nDefining field components...\n"); fflush(stdout) ;
//
//      printf("done!\nsetting up Gaussian energy...\n"); fflush(stdout) ;
//      Gaussian gauss_AB( chiAB/rho0, 2.0*a_squared, ML, A, B ) ;
//      printf("done!\nprinting to screen...\n"); fflush(stdout) ;
//      cout << "Energied! " << gauss_AB.calc_energy() << " " << U_chiab_gg << endl;
// 
//      printf("\nBefore calc all: gradwA[0][0]: %lf\n", gradwA[0][0]) ;
//
//      A.ZeroGradient() ;
//      B.ZeroGradient() ;
//      gauss_AB.calc_all() ;
//      printf("all energy: %lf force1[0][0]/rho1[0]: %lf\n", gauss_AB.energy, 
//          gauss_AB.force1[0][0]/gauss_AB.rho1[0]) ;
//      printf("gradwB: %lf force2[0][0]/rho2[0]: %lf\n", gradwB[0][0], 
//          gauss_AB.force2[0][0]/gauss_AB.rho2[0]) ;
//      printf("f1[0][0]: %lf  f2[0][0]: %lf\n", gauss_AB.force1[0][0], gauss_AB.force2[0][0]) ;
//      
//      exit(0);


      if ( myrank == 0 ) {
        printf("step %d of %d  Ubond: %lf " , step , nsteps , Ubond ) ;
//        if ( semiflex )
//          printf("Uangle: %lf ", Uangle ) ;
//        
//        if ( chiAB > 0.0 ) 
//          printf("Uchiab: %lf ", U_chiab_gg + U_chi_pg + U_chi_pp) ;
//        
//        if ( chiBC > 0.0 ) 
//          printf("Uchibc: %lf ", U_chibc);
//        
//        if ( chiAC > 0.0 ) 
//          printf("Uchiac: %lf ", U_chiac);
//
//        if ( kappa > 0.0 ) 
//          printf("Ukappa: %lf ", U_kappa_gg + U_kappa_pg + U_kappa_pp ) ;
//        
//        if ( mu != 0.0 )
//          printf("U_ms: %lf lc: %lf", U_ms, lc_order_param ) ;
        printf("\n");
        fflush( stdout ) ;
      }

      //cout << "Proc " << myrank << " now has " << ns_loc << endl; 



     // if ( step > sample_wait ) {
     //   for ( i=0 ; i<ML ; i++ )
     //     ktmp2[i] = avg_sk[0][i] / num_averages ;
 
     //   write_kspace_data( "avg_sk" , ktmp2 ) ;
     // }

  
    }// if step % print_Freq == 0
    


    ///////////////////////////////////////
    // Logarithmic spaced configurations //
    ///////////////////////////////////////
    if ( frame_freq > 0 && step == frame_freq ) {
      char nm[20] ;

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
