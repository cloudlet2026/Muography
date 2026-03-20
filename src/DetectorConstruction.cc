

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
	G4double world_sizeXY = 10. * m;
	G4double world_sizeZ = 10. * m;
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
	G4double env1_sizeX = 1.2 * m;
	G4double env1_sizeY = 1.2 * m;
	G4double env1_sizeZ = 1.2 * m;
	G4Material *env1_mat = air;
	G4double env2_sizeX = 1.2 * m;
	G4double env2_sizeY = 1.2 * m;
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


	// 圆柱体铅罐1: 高15cm, 外径10cm, 壁厚2cm, 位置(0, 0, 0)cm
	//
	G4double leadCyl1_height = 15. * cm;
	G4double leadCyl1_outerRadius = 5. * cm; // 外径10cm
	G4double leadCyl1_thickness = 2. * cm;
	G4double leadCyl1_innerRadius = leadCyl1_outerRadius - leadCyl1_thickness; // 内半径10cm
	G4Material *leadCyl_mat = plumbum;

	// 空心圆柱：直接使用G4Tubs的内外半径参数
	auto solidLeadCyl1 = new G4Tubs("LeadCyl1",
						 leadCyl1_innerRadius, leadCyl1_outerRadius, 0.5 * leadCyl1_height,
						 0. * deg, 360. * deg);

	auto logicLeadCyl1 = new G4LogicalVolume(solidLeadCyl1,
											 leadCyl_mat,
											 "LeadCyl1");
	logicLeadCyl1->SetUserLimits(userLimits);
	new G4PVPlacement(nullptr,
				  G4ThreeVector(0., 0., 0.),
				  logicLeadCyl1,
				  "LeadCyl1",
				  logicEnv1,
				  false,
				  0,
				  checkOverlaps);

	// 顶盖和底盖（实心圆柱）
	G4double leadCyl1_capThickness = leadCyl1_thickness;
	
	auto solidLeadCyl1_top = new G4Tubs("LeadCyl1_top",
						 0., leadCyl1_outerRadius, 0.5 * leadCyl1_capThickness,
						 0. * deg, 360. * deg);
	auto logicLeadCyl1_top = new G4LogicalVolume(solidLeadCyl1_top,
											 leadCyl_mat,
											 "LeadCyl1_top");
	new G4PVPlacement(nullptr,
				  G4ThreeVector(0., 0., 0.5 * leadCyl1_height + 0.5 * leadCyl1_capThickness),
				  logicLeadCyl1_top,
				  "LeadCyl1_top",
				  logicEnv1,
				  false,
				  0,
				  checkOverlaps);

	auto solidLeadCyl1_bottom = new G4Tubs("LeadCyl1_bottom",
						 0., leadCyl1_outerRadius, 0.5 * leadCyl1_capThickness,
						 0. * deg, 360. * deg);
	auto logicLeadCyl1_bottom = new G4LogicalVolume(solidLeadCyl1_bottom,
											 leadCyl_mat,
											 "LeadCyl1_bottom");
	new G4PVPlacement(nullptr,
				  G4ThreeVector(0., 0., -0.5 * leadCyl1_height - 0.5 * leadCyl1_capThickness),
				  logicLeadCyl1_bottom,
				  "LeadCyl1_bottom",
				  logicEnv1,
				  false,
				  0,
				  checkOverlaps);
	
	// 可视化属性------------------------------------------------

	auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	visAttributes->SetVisibility(false);
	logicWorld->SetVisAttributes(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(0.2, 0.4, 0.8, 0.9));
	visAttributes->SetVisibility(true);
	logicLeadCyl1->SetVisAttributes(visAttributes);
	logicLeadCyl1_top->SetVisAttributes(visAttributes);
	logicLeadCyl1_bottom->SetVisAttributes(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(0.0, 0.8, 0.5, 0.5));
	visAttributes->SetVisibility(true);
	logicEnv2->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.5));
	visAttributes->SetVisibility(true);
	logicEnv1->SetVisAttributes(visAttributes);

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
