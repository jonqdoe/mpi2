#include "pair_style.h"

#ifndef _PAIR_GAUSSERFC
#define _PAIR_GAUSSERFC


// This class is of type PairStyle and inherits all PairStyle
// routines. PairStyle is initialized first, then the Gaussian
// initialization is performed.
class GaussianErfc : public PairStyle {
  public:
    GaussianErfc() ;
    ~GaussianErfc() ;
    void Initialize_GaussianErfc(double, double, double, double, 
        int, FieldComponent, FieldComponent ) ;

};


#endif
