#include "pair_style.h"

#ifndef _PAIR_GAUSSIAN
#define _PAIR_GAUSSIAN


// This class is of type PairStyle and inherits all PairStyle
// routines. PairStyle is initialized first, then the Gaussian
// initialization is performed.
class Gaussian : public PairStyle {
  public:
    Gaussian(double, double, int, FieldComponent, FieldComponent ) ;

};


#endif
