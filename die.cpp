#include "globals.h"

void die(const char *kill) {
  fprintf(stdout,"%s\n",kill);
  fftw_mpi_cleanup();
  MPI_Finalize() ;
  exit(1);
}
