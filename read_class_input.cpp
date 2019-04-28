#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std ;

void read_class_input( const char* ) ;

int main( int argc, char** argv ) {
  read_class_input( "mpi2.input" ) ;
  return 0 ;
}





void read_class_input( const char* nm ) {

  string name(nm) ;
  ifstream inp(name) ;
  string word, line ;
  int line_no = 0 ;
  
  while ( getline(inp, line) ) {
    if ( line.length() == 0 || line.at(0) == '#')
      continue ;

    istringstream iss(line) ;
    while ( iss >> word ) {
      if ( word == "rho0" ) {
        string todouble ;
        iss >> todouble ;
        double test = stod(todouble) ;
        printf("rho0 set to be %lf\n", test ) ;
      }

      if ( word == "Nda" ) {
        string toint ;
        iss >> toint ;
        int itest = stoi( toint ) ;
        printf("Nda set to %d\n", itest ) ;
      }
    }

  }

  inp.close() ;

}
