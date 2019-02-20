#include "globals.h"

double integrate( double *inp ) {
  double l_sum = 0.0, sum = 0.0 ;
  int i ;

  for ( i=0 ; i<ML ; i++ )
    l_sum += inp[i] ;

  for ( i=0 ; i<Dim ; i++ )
    l_sum *= dx[i] ;

  MPI_Allreduce( &l_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;

  return sum ;
}
