

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"

//
#include "G4OpticalSurface.hh"
#include "G4EllipticalTube.hh"
#include "G4Tubs.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
	// 从 NIST 数据库获取常用金属材料
	// Get nist material manager
	G4NistManager *nist = G4NistManager::Instance();
	G4Material *aluminum = nist->FindOrBuildMaterial("G4_Al"); // Z=13
	G4Material *iron = nist->FindOrBuildMaterial("G4_Fe");	   // Z=26
	G4Material *tungsten = nist->FindOrBuildMaterial("G4_W");  // Z=74
	G4Material *plumbum = nist->FindOrBuildMaterial("G4_Pb");  // Z=82
	G4Material *uranium = nist->FindOrBuildMaterial("G4_U");   // Z=92
	G4Material *air = nist->FindOrBuildMaterial("G4_AIR");

	// Option to switch on/off checking of volumes overlaps
	//
	G4bool checkOverlaps = false;

	//
	// World
	//
	G4double world_sizeXY = 1.2 * m;
	G4double world_sizeZ = 2. * m;
	G4Material *world_mat = air;

	auto solidWorld = new G4Box("World",													// its name
								0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ); // its size

	auto logicWorld = new G4LogicalVolume(solidWorld, // its solid
										  world_mat,  // its material
										  "World");	  // its name

	// 控制最大步长，避免跨界步过大导致命中漂移
	auto userLimits = new G4UserLimits(0.1 * mm);
	logicWorld->SetUserLimits(userLimits);

	auto physWorld = new G4PVPlacement(nullptr,			// no rotation
									   G4ThreeVector(), // at (0,0,0)
									   logicWorld,		// its logical volume
									   "World",			// its name
									   nullptr,			// its mother  volume
									   false,			// no boolean operation
									   0,				// copy number
									   checkOverlaps);	// overlaps checking
	G4double env1_sizeX = world_sizeXY;
	G4double env1_sizeY = world_sizeXY;
	G4double env1_sizeZ = 1.2 * m;
	G4Material *env1_mat = air;
	G4double env2_sizeX = world_sizeXY;
	G4double env2_sizeY = world_sizeXY;
	G4double env2_sizeZ = 1.6 * m;
	G4Material *env2_mat = air;
	//
	// Envelope2
	//
	auto solidEnv2 = new G4Box("Envelope2",											  // its name
							   0.5 * env2_sizeX, 0.5 * env2_sizeY, 0.5 * env2_sizeZ); // its size

	auto logicEnv2 = new G4LogicalVolume(solidEnv2,	   // its solid
										 env2_mat,	   // its material
										 "Envelope2"); // its name
	logicEnv2->SetUserLimits(userLimits);
	new G4PVPlacement(nullptr,		   // no rotation
					  G4ThreeVector(), // at (0,0,0)
					  logicEnv2,	   // its logical volume
					  "Envelope2",	   // its name
					  logicWorld,	   // its mother  volume
					  false,		   // no boolean operation
					  0,			   // copy number
					  checkOverlaps);  // overlaps checking
	//
	// Envelope1
	//
	auto solidEnv1 = new G4Box("Envelope1",											  // its name
							   0.5 * env1_sizeX, 0.5 * env1_sizeY, 0.5 * env1_sizeZ); // its size

	auto logicEnv1 = new G4LogicalVolume(solidEnv1,	   // its solid
										 env1_mat,	   // its material
										 "Envelope1"); // its name
	logicEnv1->SetUserLimits(userLimits);
	new G4PVPlacement(nullptr,		   // no rotation
					  G4ThreeVector(), // at (0,0,0)
					  logicEnv1,	   // its logical volume
					  "Envelope1",	   // its name
					  logicEnv2,	   // its mother  volume
					  false,		   // no boolean operation
					  0,			   // copy number
					  checkOverlaps);  // overlaps checking


	//
	// 空心铅盒 (2cm厚, 外部尺寸10x10x10 cm³), 放在铁箱子内部
	//
	G4double leadBox_outerSize = 10. * cm;   // 外部尺寸 10cm
	G4double leadBox_thickness = 2. * cm;     // 壁厚 2cm
	G4double leadBox_innerSize = leadBox_outerSize - 2. * leadBox_thickness;  // 内部尺寸 6cm
	G4Material *leadBox_mat = plumbum;

	// 外壳
	auto solidLeadBoxOuter = new G4Box("LeadBoxOuter",
										0.5 * leadBox_outerSize, 0.5 * leadBox_outerSize, 0.5 * leadBox_outerSize);
	// 内腔
	auto solidLeadBoxInner = new G4Box("LeadBoxInner",
										0.5 * leadBox_innerSize, 0.5 * leadBox_innerSize, 0.5 * leadBox_innerSize);
	// 空心盒
	auto solidLeadBox = new G4SubtractionSolid("LeadBox", solidLeadBoxOuter, solidLeadBoxInner);

	auto logicLeadBox = new G4LogicalVolume(solidLeadBox,
											leadBox_mat,
											"LeadBox");
	logicLeadBox->SetUserLimits(userLimits);
	new G4PVPlacement(nullptr,
					  G4ThreeVector(0, 0, 0),  // 放在铁箱子中心（即 Envelope1 中心）
					  logicLeadBox,
					  "LeadBox",
					  logicEnv1,  // 母体是 Envelope1（和铁箱子共享同一母体）
					  false,
					  0,
					  checkOverlaps);
	// 可视化属性------------------------------------------------

	auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	visAttributes->SetVisibility(false);
	logicWorld->SetVisAttributes(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(0.0, 0.8, 0.5, 0.5));
	visAttributes->SetVisibility(true);
	logicEnv2->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.5));
	visAttributes->SetVisibility(true);
	logicEnv1->SetVisAttributes(visAttributes);

	// 铅盒可视化
	visAttributes = new G4VisAttributes(G4Colour(0.7, 0.5, 0.9, 0.5));
	visAttributes->SetVisibility(true);
	logicLeadBox->SetVisAttributes(visAttributes);

	if (logicShape1)
	{
		visAttributes = new G4VisAttributes(G4Colour(0.0, 0.66, 0., 0.5));
		visAttributes->SetVisibility(true);
		logicShape1->SetVisAttributes(visAttributes);
	}
	if (logicShape2)
	{
		visAttributes = new G4VisAttributes(G4Colour(0.0, 0.66, 0., 0.5));
		visAttributes->SetVisibility(true);
		logicShape2->SetVisAttributes(visAttributes);
	}
	if (logicShape3)
	{
		visAttributes = new G4VisAttributes(G4Colour(0.0, 0.66, 0., 0.5));
		visAttributes->SetVisibility(true);
		logicShape3->SetVisAttributes(visAttributes);
	}
	if (logicShape4)
	{
		visAttributes = new G4VisAttributes(G4Colour(0.0, 0.66, 0., 0.5));
		visAttributes->SetVisibility(true);
		logicShape4->SetVisAttributes(visAttributes);
	}
	//
	// always return the physical World
	//
	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
