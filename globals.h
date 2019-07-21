#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "field_component.h"
#include "pair_style_gaussian.h"
#include "pair_style_gauss-erfc.h"
#include "pair_style_erfc2.h"
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
FieldComponent *Components ;

#ifndef MAIN
extern
#endif
Gaussian *Gauss ;

#ifndef MAIN
extern
#endif
GaussianErfc *GaussErfc ;

#ifndef MAIN
extern
#endif
ErfcErfc *Erfc2 ;


#ifndef MAIN
extern
#endif
double **x, **xtmp, **f , *tmp, **grid_W, V, L[Dim], dx[Dim], gvol, Lh[Dim], 
       **rho , **w ,  *r_dudr , *tmp2, 
       *Diff, delt , Pvir, Ptens[Dim][Dim], Unb, *mass,
       Stress_bonds[Dim][Dim], mem_use,
       *anneal_chi,
       z_min, z_max, send_buff, **rec_N_x, **rec_S_x, **send_S_x, **send_N_x,
       *bond_coeff, *bond_eq, Ubond,
       **angle_coeff, Uangle,
       *gaussian_prefactor, *gaussian_sigma,
       *gausserfc_prefactor, *gausserfc_sigma, *gausserfc_Rp,
       *gausserfc_xi,
       *erfc2_prefactor, *erfc2_Rp, *erfc2_xi ;

#ifndef MAIN
extern
#endif
int nstot, *molecID, *tp, 
    Nx[Dim], M, nsteps, step, print_freq , 
    **grid_inds, pmeorder, spline_weights, lagrange_weights, 
    grid_per_partic, ntypes , stress_freq , buff_size, buff_ind ,
    sample_wait, sample_freq,
    frame_freq, traj_freq, do_anneal, *anneal_steps, n_anneal_pts, next_anneal_update,
    cur_anneal_step,
    NxL[Dim], ML, zstart, size, nprocs, myrank,
    *my_inds, ns_loc, n_N_ghost, n_S_ghost, *send_N_inds, *send_S_inds, 
    Sbound_ct, Nbound_ct, *rec_N_inds, *rec_S_inds, total_ghost, *ghost_inds,
    *n_bonds, **bonded_to, **bond_type, n_total_bonds,
    *n_angles, **angle_first, **angle_mid, **angle_end, **angle_type, n_total_angles, init_flag,
    *local_flag,
    nbond_types, nangle_types,
    n_gaussian_pairstyles, **gaussian_types,
    n_gausserfc_pairstyles, **gausserfc_types,
    n_erfc2_pairstyles, **erfc2_types ;


#ifndef MAIN
extern
#endif
char **xc, tt[80] ;


#ifndef MAIN
extern
#endif
complex<double> *ktmp2, *ktmp, I, **avg_sk ;

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
