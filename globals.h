#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "mpi.h"
#include "fftw3-mpi.h"
#include <Eigen/Dense>

using namespace std ;
using namespace Eigen ;

#define PI   3.141592653589793238462643383

#define KDelta(i,j) ( i==j ? 1.0 : 0.0 ) 
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define min(A,B) ((A)<(B) ? (A) : (B) )

#define Dim 2


#ifndef MAIN
extern
#endif
double **x, **xtmp, **f , *tmp, **grid_W, V, L[Dim], dx[Dim], gvol, Lh[Dim], 
       **rho , **w ,  *r_dudr , *tmp2, 
       *Diff, delt , chiAC, chiBC, chiAB, kappa, rho0 , num_averages , 
       Ubond , Uchi, Ukappa, Ptens[Dim][Dim], Pvir , phiHA, phiHB, phiHC,
       *rhoha, *rhohb, *rhohc, *rhoda, *rhodb , *rhot, *rhop, *smrhop,
       **rhoha_t, **rhohc_t, **rhohb_t, **rhoda_t, **rhodb_t , **rhop_t,
       Stress_bonds[Dim][Dim],
       ***sts_buf , Rg, Rg3, phiP, Rp, Xi, Vp, *gammaP,
       *gradwA[Dim], *gradwB[Dim], *gradwP[Dim], *gradwC[Dim],
       *uG, *grad_uG[Dim], *uP, *grad_uP[Dim], *uPG, *grad_uPG[Dim] , mem_use,
       U_chiab_gg, U_chi_pg, U_chi_pp, U_kappa_gg, U_kappa_pg, U_kappa_pp,
       U_chibc, U_chiac,
       *phiA, *Dweights, 
       *anneal_chi, a_squared ,
       z_min, z_max, send_buff, **rec_N_x, **rec_S_x, **send_S_x, **send_N_x,
       **bond_coeff, **bond_eq, 
       **angle_coeff, semiflex_k, semiflex_req, semiflex_lam, Uangle,
       **mono_u, ***mono_S, mu, ***S_field, ***S_conv_u, ***tmp_tensor, U_ms,
       lc_order_param,
       *gamma_sig, *uAG, *grad_uAG[Dim], eps ; // Parameters for attractive nanoparticles

#ifndef MAIN
extern
#endif
int nstot, *tp, nC, Nhc, nA, nB, nD, Nha, Nhb, Nda, Ndb, nP, 
    Nx[Dim], M, nsteps, step, print_freq , nsC, nsD, nsA, nsB, 
    **grid_inds, pmeorder, spline_weights, lagrange_weights, 
    grid_per_partic, ntypes , stress_freq , buff_size, buff_ind ,
    sample_wait, sample_freq, B_partics , 
    nthreads, frame_freq, traj_freq, do_anneal, *anneal_steps, n_anneal_pts, next_anneal_update,
    cur_anneal_step,
    NxL[Dim], ML, zstart, size, nprocs, myrank,
    *my_inds, ns_loc, n_N_ghost, n_S_ghost, *send_N_inds, *send_S_inds, 
    Sbound_ct, Nbound_ct, *rec_N_inds, *rec_S_inds, total_ghost, *ghost_inds,
    *n_bonds, **bonded_to,
    semiflex, *n_angles, **angle_first, **angle_mid, **angle_end, n_tot_angles, init_flag,
    *local_flag ;


#ifndef MAIN
extern
#endif
char **xc, tt[80] ;


#ifndef MAIN
extern
#endif
complex<double> *ktmp2, *ktmp, I, **avg_sk , *grad_uG_hat[Dim] , *grad_uP_hat[Dim], *grad_uPG_hat[Dim],
  *grad_uAG_hat[Dim] ;

#ifndef MAIN
extern
#endif
long idum;



#ifndef MAIN
extern
#endif
fftw_complex *fmin0, *fmot0 ;

#ifndef MAIN
extern
#endif
fftw_plan fwd0, fbk0 ;







double ran2(void ) ;
double gasdev2( void ) ;
int cell_stack( int* );
void cell_unstack( int , int* );
void field_gradient( double* , double* , int ) ;
void field_gradient_cdif( double* , double* , int ) ;
void convolve_fields( double*, double*, double* ) ;
void write_grid_data( const char* , double* ) ;
void write_Sfield_data( const char* , double*** ) ;
void write_kspace_data( const char* , complex<double>* ) ;

int unstack_stack( int ) ;
void unstack_local( int, int* ) ;
int stack( int* ) ;
int stack_to_local( int* ) ;
int stack_local( int* ) ;
double integrate( double* );
void unstack( int , int* );
double get_k( int , double* ) ;
double get_k_alias( int , double* ) ;
void get_r( int , double* ) ;

void fftw_fwd( double* , complex<double>* );
void fftw_back( complex<double>* , double* );
int remove_dupes( int* , int );

int malloc2ddouble( double***, int , int  ) ;
int malloc2dint( int***, int , int  ) ;

void pbc_vdr( double*, double* , double* );
double pbc_mdr2( double*, double* , double* );
void die(const char *kill);

int mpi_Bcast_posits( int, int, double**, int, MPI_Comm ) ;
int allocate_mpi_tmp( int, int ) ;
void find_passengers( void ) ;
int is_partic_in_list( int, int, int* ) ;

void add_to_array( int, int*, double**, double** ) ;
void insert_in_array( int, int*, double**, double** ) ;
void extract_from_array( int, int*, double**, double** ) ;

double double_dot( double**, double** ) ;
int particle_has_orientation( int ) ;
double  diag_mat( double**, double* ) ;
