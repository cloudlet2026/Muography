

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Cache.hh"
#include "globals.hh"
#include <CLHEP/Units/SystemOfUnits.h>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4VisAttributes;
class G4Material;

class Materials;
class PhotonDetSD;
// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction() = default;
  ~DetectorConstruction() override;

  G4VPhysicalVolume *Construct() override;

private:
  G4LogicalVolume *logicShape1 = nullptr;
  G4LogicalVolume *logicShape2 = nullptr;
  G4LogicalVolume *logicShape3 = nullptr;
  G4LogicalVolume *logicShape4 = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
