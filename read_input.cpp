#include "globals.h"
#include <sstream>
#include <string>

void read_anneal( void ) ;
void set_defaults( void ) ;
void write_runtime_parameters( void ) ;


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
      if ( word == "Nda" ) {
        string toint ;
        iss >> toint ;
        Nda = stoi( toint ) ;
      }
      else if ( word == "Ndb" ) {
        string toint ;
        iss >> toint ;
        Ndb = stoi( toint ) ;
      }

      else if ( word == "phiha" ) {
        string todouble ;
        iss >> todouble ;
        phiHA = stod(todouble) ;
      }
      else if ( word == "phihb" ) {
        string todouble ;
        iss >> todouble ;
        phiHB = stod(todouble) ;
      }
      else if ( word == "phihc" ) {
        string todouble ;
        iss >> todouble ;
        phiHC = stod(todouble) ;
      }

      else if ( word == "Nha" ) {
        string toint ;
        iss >> toint ;
        Nha = stoi( toint ) ;
      }
      else if ( word == "Nhb" ) {
        string toint ;
        iss >> toint ;
        Nhb = stoi( toint ) ;
      }
      else if ( word == "Nhc" ) {
        string toint ;
        iss >> toint ;
        Nhc = stoi( toint ) ;
      }

      else if ( word == "smear_a" ) {
        string todouble ;
        iss >> todouble ;
        a_squared = stod(todouble) ;
        a_squared *= a_squared ;
      }

      else if ( word == "semiflex" ) {
        string toint ;
        iss >> toint ;
        semiflex = stoi( toint ) ;
      }
      else if ( word == "semiflex_req" ) {
        string todouble ;
        iss >> todouble ;
        semiflex_req = stod(todouble) ;
      }
      else if ( word == "semiflex_k" ) {
        string todouble ;
        iss >> todouble ;
        semiflex_k = stod(todouble) ;
      }
      else if ( word == "semiflex_lam" ) {
        string todouble ;
        iss >> todouble ;
        semiflex_lam = stod(todouble) ;
      }

      else if ( word == "phiP" ) {
        string todouble ;
        iss >> todouble ;
        phiP = stod(todouble) ;
      }
      else if ( word == "Rp" ) {
        string todouble ;
        iss >> todouble ;
        Rp = stod(todouble) ;
      }
      else if ( word == "Xi" ) {
        string todouble ;
        iss >> todouble ;
        Xi = stod(todouble) ;
      }
      else if ( word == "B_partics" ) {
        string toint ;
        iss >> toint ;
        B_partics = stoi( toint ) ;
      }
      else if ( word == "eps" ) {
        string todouble ;
        iss >> todouble ;
        eps = stod(todouble) ;
      }

      else if ( word == "rho0" ) {
        string todouble ;
        iss >> todouble ;
        rho0 = stod(todouble) ;
      }
      else if ( word == "chiAB" ) {
        string todouble ;
        iss >> todouble ;
        chiAB = stod(todouble) ;
      }
      else if ( word == "chiAC" ) {
        string todouble ;
        iss >> todouble ;
        chiAC = stod(todouble) ;
      }
      else if ( word == "chiBC" ) {
        string todouble ;
        iss >> todouble ;
        chiBC = stod(todouble) ;
      }
      else if ( word == "kappa" ) {
        string todouble ;
        iss >> todouble ;
        kappa = stod(todouble) ;
      }

      else if ( word == "DiffA" ) {
        string todouble ;
        iss >> todouble ;
        Diff[0] = stod(todouble) ;
      }
      else if ( word == "DiffB" ) {
        string todouble ;
        iss >> todouble ;
        Diff[1] = stod(todouble) ;
      }
      else if ( word == "DiffC" ) {
        string todouble ;
        iss >> todouble ;
        Diff[2] = stod(todouble) ;
      }

      else if ( word == "Lx" ) {
        string todouble ;
        iss >> todouble ;
        L[0] = stod(todouble) ;
      }
      else if ( word == "Ly" ) {
        string todouble ;
        iss >> todouble ;
        L[1] = stod(todouble) ;
      }
      else if ( word == "Lz" ) {
        if ( Dim == 2 ) {
          cout << "Lz parameter ignored for 2D simulation!\n" ;
        }
        else {
          string todouble ;
          iss >> todouble ;
          L[2] = stod(todouble) ;
        }
      }

      else if ( word == "Nx" ) {
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

    }//line stream

  }//getline()

  inp.close() ;

  write_runtime_parameters() ;

  read_anneal() ;

}

void set_defaults() {

  Nda = 20 ;
  Ndb = 20 ;
  phiHA = phiHB = phiHC = 0.0 ;
  Nha = Nhb = Nhc = 5 ;
  a_squared = 1.0 ;

  semiflex = 0 ;
  semiflex_req = 1.0 ;
  semiflex_k = 10.0 ;
  semiflex_lam = 1.0 ;

  phiP = 0.0 ;
  Rp = 1.0 ;
  Xi = 0.1 ;
  B_partics = 0 ;
  eps = 0.0 ;

  rho0 = 5.0 ;
  chiAB = 0.5 ;
  chiAC = chiBC = 0.0 ;

  kappa = 10.0 ;

  Diff[0] = Diff[1] = Diff[2] = 1.0 ;
  
  for (int j=0 ; j<Dim ; j++ ) {
    L[j] = 20.0 ;
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

  otp << "Nda " << Nda << endl;
  otp << "Ndb " << Ndb << endl;
  otp << "phiha " << phiHA << endl;
  otp << "phihb " << phiHB << endl;
  otp << "phihc " << phiHC << endl;
  otp << "Nha " << Nha << endl;
  otp << "Nhb " << Nhb << endl;
  otp << "Nhc " << Nhc << endl;
  otp << "smear_a " << sqrt(a_squared) << endl;
  otp << "semiflex " << semiflex << endl;
  otp << "semiflex_req " << semiflex_req << endl;
  otp << "semiflex_k " << semiflex_k << endl;
  otp << "semiflex_lam " << semiflex_lam << endl;
  otp << "phiP " << phiP << endl;
  otp << "Rp " << Rp << endl;
  otp << "Xi " << Xi << endl;
  otp << "B_partics " << B_partics << endl;
  otp << "eps " << eps << endl;
  otp << "rho0 " << rho0 << endl;
  otp << "chiAB " << chiAB << endl;
  otp << "chiAC " << chiAC << endl;
  otp << "chiBC " << chiBC << endl;
  otp << "kappa " << kappa << endl;
  otp << "DiffA " << Diff[0] << endl;
  otp << "DiffB " << Diff[1] << endl;
  otp << "DiffC " << Diff[2] << endl;
  otp << "Lx " << L[0] << endl;
  otp << "Ly " << L[1] << endl;
  otp << "Nx " << Nx[0] << endl;
  otp << "Ny " << Nx[1] << endl;
  if ( Dim == 3 ) {
    otp << "Lz " << L[2] << endl;
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

  otp.close() ;
}
