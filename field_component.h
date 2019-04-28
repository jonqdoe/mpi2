#ifndef _FIELD_COMP
#define _FIELD_COMP

class FieldComponent {
  public:
    double *rho, **gradU ;

    void ZeroGradient() ;
    FieldComponent(int) ;
    ~FieldComponent() ;
} ;

#endif
