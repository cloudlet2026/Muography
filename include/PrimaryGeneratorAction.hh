#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <TF2.h>
#include <TFile.h>
#include <TH1.h>

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  // method from the base class
  void GeneratePrimaries(G4Event *) override;

  // method to access particle gun
  const G4GeneralParticleSource *GetParticleGun() const { return fParticleGun; }

private:
  // G4ParticleGun *fParticleGun = nullptr; // pointer a to G4 gun class
  G4GeneralParticleSource *fParticleGun = nullptr;
  G4Box *fEnvelopeBox = nullptr;
  static G4bool fSeaLevel;
  TF2 *fFlux = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
