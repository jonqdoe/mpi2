
class Molecule {
  public:
    double *rho ; // Density field for this molecule
    int Ns ; // Sites per molecule
    int nmolecs ; // Number of molecules of this type

    Molecule() ;
    ~Molecule() ;

} ;
