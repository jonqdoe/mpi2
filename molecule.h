
class Molecule {
  public:
    int *field_type ; // Component field for this molecule type
    int Ns ; // Sites per molecule
    int nmolecs ; // Number of molecules of this type

    Molecule(int) ;
    ~Molecule() ;

} ;
