#ifndef G4SBSGlobalField_hh
#define G4SBSGlobalField_hh

#include "globals.hh"
#include "sbstypes.hh"
#include <vector>

#include "G4SBSMagneticField.hh"

class TH2F;

class G4SBSGlobalField : public G4MagneticField {
public:
  G4SBSGlobalField();
  ~G4SBSGlobalField();

  void GetFieldValue( const  double Point[3], double *Bfield ) const;

  void SetInvertField( G4bool b );

  void AddToscaField(const char *); //return a pointer to the added field so that we can store it in G4SBSDetectorConstruction

  void AddField( G4SBSMagneticField *f );
  void DropField( G4SBSMagneticField *f );

  void DebugField();

  void ScaleFields( G4double, Arm_t );
  
  std::vector<TH2F *> fFieldPlots;

  bool fInverted;
  
  private:
  std::vector<G4SBSMagneticField *> fFields;
  //std::vector<G4double> fScaleFactor;
};

#endif//G4SBSGlobalField_hh
