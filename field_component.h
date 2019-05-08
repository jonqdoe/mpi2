#ifndef _FIELD_COMP
#define _FIELD_COMP

class FieldComponent {
  public:
    double *rho, **force ;

    void ZeroGradient() ;
    void Initialize(int) ;
    FieldComponent() ;
    ~FieldComponent() ;
} ;

#endif
