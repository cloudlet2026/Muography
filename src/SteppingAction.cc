#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Run.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UImanager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UnitsTable.hh"

// Purpose: Save relevant information into User Track Information

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction *detector,
							   EventAction *event)
	: fDetector(detector), fEventAction(event)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *theStep)
{
	if (theStep->GetTrack()->GetParticleDefinition()->GetParticleName() != "mu-")
	{
		return;
	}

	G4StepPoint *thePrePoint = theStep->GetPreStepPoint();
	G4StepPoint *thePostPoint = theStep->GetPostStepPoint();

	// 只在真正跨越几何边界时记录
	if (thePostPoint->GetStepStatus() != fGeomBoundary)
	{
		return;
	}

	G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
	G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();

	G4String thePrePVname = " ";
	G4String thePostPVname = " ";

	if (!thePrePV || !thePostPV)
	{
		return;
	}
	thePrePVname = theStep->GetTrack()->GetVolume()->GetName();
	thePostPVname = theStep->GetTrack()->GetNextVolume()->GetName();
	// thePrePVname = thePrePV->GetName();
	// thePostPVname = thePostPV->GetName();

	const G4ThreeVector &prePos = thePrePoint->GetPosition();
	const G4ThreeVector &postPos = thePostPoint->GetPosition();
	G4ThreeVector dir = postPos - prePos;

	auto projectToZ = [&](G4double z_target) -> G4ThreeVector
	{
		if (std::abs(dir.z()) < 1e-12)
		{
			return postPos; // 方向几乎水平，直接返回post位置
		}
		double t = (z_target - prePos.z()) / dir.z();
		return prePos + t * dir;
	};

	// 实际几何：World -> Envelope2 -> Envelope1 -> Shape1
	const G4double z_env2 = 0.8 * m;
	const G4double z_env1 = 0.6 * m;

	if (thePrePVname == "World" && thePostPVname == "Envelope2")
	{
		fEventAction->AddBoundaryCrossing(1, projectToZ(z_env2));
	}
	else if (thePrePVname == "Envelope2" && thePostPVname == "Envelope1")
	{
		fEventAction->AddBoundaryCrossing(2, projectToZ(z_env1));
	}
	else if (thePrePVname == "Envelope1" && thePostPVname == "Envelope2")
	{
		fEventAction->AddBoundaryCrossing(3, projectToZ(-z_env1));
	}
	else if (thePrePVname == "Envelope2" && thePostPVname == "World")
	{
		fEventAction->AddBoundaryCrossing(4, projectToZ(-z_env2));
	}
}
