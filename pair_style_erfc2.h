#include "pair_style.h"

#ifndef _PAIR_ERFCERFC
#define _PAIR_ERFCERFC


// This class is of type PairStyle and inherits all PairStyle
// routines. PairStyle is initialized first, then the Gaussian
// initialization is performed.
class ErfcErfc : public PairStyle {
  public:
    ErfcErfc() ;
    ~ErfcErfc() ;
    void Initialize_Erfc2(double, double, double, 
        int, FieldComponent, FieldComponent ) ;

};


#endif
