#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4GenericMessenger.hh"
#include "G4OpticalSurface.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalSkinSurface.hh"
#include "Randomize.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"




#include "detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	MyDetectorConstruction();
	~MyDetectorConstruction();
	
	virtual G4VPhysicalVolume *Construct();
	
	G4ThreeVector GenerateVertex(const G4String& region) const;
	
private:
	G4LogicalVolume *logicDetector,*photocathode_logic;
	virtual void ConstructSDandField();
	G4OpticalSurface* GetPhotOptSurf();
	
	G4double DetPos,zGun,transparency;
	G4bool CaseA,CaseB,CaseC,CaseD,PMTR,Gas,Patron1,Patron2,Patron3;
	
	G4Box *solidWorld;
	G4Tubs *solidTub1,*solidTub2,*solidDetector;
	G4Tubs *solidTub2Teflon,*solidTub1TeflonInt,*solidTub1TeflonExt,*solidMesh;
	G4Tubs *solidDetectorTapa,*solidGasfilm;
	
	G4Sphere *solidSphere,*solidSphereTeflon;
	G4LogicalVolume *logicWorld,*logicSphere,*logicSphereTeflon,*logicTub1,*logicTub2,*logicTub2Teflon;
	G4VPhysicalVolume *physWorld,*physSphere,*physSphereTeflon,*physTub1,*physTub2,*physDetector,*physTub2Teflon;
	
	G4LogicalVolume *logicTub1TeflonInt,*logicTub1TeflonExt,*logicMesh,*logicDetectorTapa,*logicGasfilm;
	G4VPhysicalVolume *physTub1TeflonInt,*physTub1TeflonExt,*physMesh,*physDetectorTapa,*physGasfilm;
	
	G4Material *lXe,*Steel,*Teflon,*Kovar,*silica,*FakeGrid,*GXe;
	
	
	G4GenericMessenger *fMessenger;
	
	void DefineMaterials();
	
};

#endif
