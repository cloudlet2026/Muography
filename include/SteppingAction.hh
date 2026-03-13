#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;
class G4LogicalVolume;

class G4OpBoundaryProcess;
class G4Track;
class G4StepPoint;

/// Stepping action class

class SteppingAction : public G4UserSteppingAction
{
public:
	SteppingAction(DetectorConstruction *, EventAction *);
	~SteppingAction() override = default;

	// method from the base class
	void UserSteppingAction(const G4Step *) override;

private:
	G4int EventID = 0;
	DetectorConstruction *fDetector = nullptr;
	EventAction *fEventAction = nullptr;

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
};
#endif
