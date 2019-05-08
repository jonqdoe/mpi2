#include "globals.h"
#include <sstream>
#include <string>

void read_anneal( void ) ;
void set_defaults( void ) ;
void write_runtime_parameters( void ) ;
void allocate_gaussians( void ) ;

void read_input( void ) {

  int i,j;
  double d1 ;
  char tt[80] ;
  Diff = ( double* ) calloc( 4 , sizeof( double ) ) ;

  set_defaults() ;

  ifstream inp("mpi2.input") ;
  string word, line ;

  while ( getline(inp, line) ) {
    // Blank or commented line
    if ( line.length() == 0 || line.at(0) == '#' )
      continue ;

    istringstream iss(line) ;
    // Loop over words in line
    while ( iss >> word ) {

      if ( word == "Nx" ) {
        string toint ;
        iss >> toint ;
        Nx[0] = stoi( toint ) ;
      }
      else if ( word == "Ny" ) {
        string toint ;
        iss >> toint ;
        Nx[1] = stoi( toint ) ;
      }
      else if ( word == "Nz" ) {
        if ( Dim == 2 ) {
          cout << "Nz parameter ignored for 2D simulation!\n" ;
        }
        else {
          string toint ;
          iss >> toint ;
          Nx[2] = stoi( toint ) ;
        }
      }

      else if ( word == "delt" ) {
        string todouble ;
        iss >> todouble ;
        delt = stod(todouble) ;
      }
      else if ( word == "pmeorder" ) {
        string toint ;
        iss >> toint ;
        pmeorder = stoi( toint ) ;
      }
      else if ( word == "send_buff" ) {
        string todouble ;
        iss >> todouble ;
        send_buff = stod(todouble) ;
      }

      else if ( word == "nsteps" ) {
        string toint ;
        iss >> toint ;
        nsteps = stoi( toint ) ;
      }
      else if ( word == "print_freq" ) {
        string toint ;
        iss >> toint ;
        print_freq = stoi( toint ) ;
      }
      else if ( word == "sample_wait" ) {
        string toint ;
        iss >> toint ;
        sample_wait = stoi( toint ) ;
      }
      else if ( word == "sample_freq" ) {
        string toint ;
        iss >> toint ;
        sample_freq = stoi( toint ) ;
      }
      else if ( word == "stress_freq" ) {
        string toint ;
        iss >> toint ;
        stress_freq = stoi( toint ) ;
      }
      else if ( word == "frame_freq" ) {
        string toint ;
        iss >> toint ;
        frame_freq = stoi( toint ) ;
      }
      else if ( word == "traj_freq" ) {
        string toint ;
        iss >> toint ;
        traj_freq = stoi( toint ) ;
      }
      else if ( word == "init_flag" ) {
        string toint ;
        iss >> toint ;
        init_flag = stoi( toint ) ;
      }

      // Gaussian interactions //
      else if ( word == "n_gaussians" ) {
        string toint ;
        iss >> toint ;
        n_gaussian_pairstyles = stoi( toint ) ;
        allocate_gaussians() ;

        for ( i=0 ; i<n_gaussian_pairstyles; i++ ) {
          getline( inp, line ) ;
          istringstream iss2(line) ;

          string convert ;
          iss2 >> convert ;
          if ( convert != "gaussian" ) {
            die("Error in Gaussian stuff in input file!\n");
          }

          iss2 >> convert ;
          gaussian_types[i][0] = stoi( convert ) ;
          
          iss2 >> convert ;
          gaussian_types[i][1] = stoi( convert ) ;
          
          iss2 >> convert ;
          gaussian_prefactor[i] = stod( convert ) ;

          iss2 >> convert ;
          gaussian_sigma[i] = stod( convert ) ;
        }
      }

    }//line stream

  }//getline()

  inp.close() ;

  write_runtime_parameters() ;

  read_anneal() ;

}

void set_defaults() {

  for (int j=0 ; j<Dim ; j++ ) {
    Nx[j] = 35 ;
  }

  delt = 0.002 ;
  pmeorder = 2 ;
  send_buff = 2.5 ;

  nsteps = 250000 ;
  print_freq = 500 ;
  sample_wait = nsteps + 1 ;
  sample_freq = 5 ;

  stress_freq = 0 ;
  frame_freq = 0 ;

  traj_freq = 1000 ;
  init_flag = 0 ;
}


void write_runtime_parameters() {
  ofstream otp ;
  otp.open("tild_out.input" ) ;

  otp << "Nx " << Nx[0] << endl;
  otp << "Ny " << Nx[1] << endl;
  if ( Dim == 3 ) {
    otp << "Nz " << endl;
  }
  otp << "delt " << delt << endl;
  otp << "pmeorder " << pmeorder << endl;
  otp << "send_buff " << send_buff << endl;
  otp << "nsteps " << nsteps << endl;
  otp << "print_freq " << print_freq << endl;
  otp << "sample_wait " << sample_wait << endl;
  otp << "sample_freq " << sample_freq << endl;
  otp << "stress_freq " << stress_freq << endl;
  otp << "frame_freq " << frame_freq << endl;
  otp << "traj_freq " << traj_freq << endl;
  otp << "init_flag " << init_flag << endl;

  otp << "n_gaussians" << n_gaussian_pairstyles << endl;
  for ( int i=0 ; i<n_gaussian_pairstyles ; i++ ) {
    otp << "gaussian " << gaussian_types[i][0] << " " << gaussian_types[i][1] \
      << " " << gaussian_prefactor[i] << " " << gaussian_sigma[i] << endl; 
  }

  otp.close() ;
}
