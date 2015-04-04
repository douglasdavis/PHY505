//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fCheckOverlaps(true),
    fNofLayers(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_WATER");
  nistManager->FindOrBuildMaterial("G4_POLYPROPYLENE");
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
		 kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers             = 1;
  G4double layThickness  = 10.*cm;
  G4double lay1Thickness = 10.*cm;
  G4double lay2Thickness = 10.*cm;
  G4double lay3Thickness = 10.*cm;
  G4double lay4Thickness = 10.*cm;
  G4double lay5Thickness = 10.*cm;
  G4double lay6Thickness = 10.*cm;
  G4double lay7Thickness = 10.*cm;
  G4double lay8Thickness = 10.*cm;
  
  G4double calorSizeXY   = 10.*cm;

  G4double layersThickness = 8*layThickness;

  G4double calorThickness = fNofLayers * layersThickness;
  G4double worldSizeXY    = 1.5 * calorSizeXY;
  G4double worldSizeZ     = 1.5 * calorThickness; 
  
  // Get materials
  std::string laymat_Pb  = "G4_Pb";
  std::string laymat_PP  = "G4_POLYPROPYLENE";
  std::string laymat_H2O = "G4_WATER";
  
  std::string laymat = laymat_H2O;
  G4Material* defaultMaterial   = G4Material::GetMaterial("Galactic");
  G4Material* lay1Material      = G4Material::GetMaterial(laymat);
  G4Material* lay2Material      = G4Material::GetMaterial(laymat);
  G4Material* lay3Material      = G4Material::GetMaterial(laymat);
  G4Material* lay4Material      = G4Material::GetMaterial(laymat);
  G4Material* lay5Material      = G4Material::GetMaterial(laymat);
  G4Material* lay6Material      = G4Material::GetMaterial(laymat);
  G4Material* lay7Material      = G4Material::GetMaterial(laymat);
  G4Material* lay8Material      = G4Material::GetMaterial(laymat);

  
  if ( ! defaultMaterial || ! lay1Material || ! lay2Material ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()",
		"MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  G4VSolid* worldS = new G4Box("World",worldSizeXY/2,worldSizeXY/2,worldSizeZ/2); // its size
                         
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS,           // its solid
						 defaultMaterial,  // its material
						 "World");         // its name
                                   
  G4VPhysicalVolume* worldPV = new G4PVPlacement(0,                // no rotation
						 G4ThreeVector(),  // at (0,0,0)
						 worldLV,          // its logical volume                         
						 "World",          // its name
						 0,                // its mother  volume
						 false,            // no boolean operation
						 0,                // copy number
						 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  G4VSolid* calorimeterS = new G4Box("Calorimeter",     // its name
				     calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                         
  G4LogicalVolume* calorLV = new G4LogicalVolume(calorimeterS,     // its solid
						 defaultMaterial,  // its material
						 "Calorimeter");   // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(),  // at (0,0,0)
		    calorLV,          // its logical volume                         
		    "Calorimeter",    // its name
		    worldLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 
   
  //                                 
  // Layer
  //
  
  G4VSolid* layerS = new G4Box("Layer",           // its name
			       calorSizeXY/2, calorSizeXY/2, layersThickness/2); //its size
                         
  G4LogicalVolume* layerLV = new G4LogicalVolume(layerS,           // its solid
						 defaultMaterial,  // its material
						 "Layer");         // its name

  new G4PVReplica("Layer",          // its name
		  layerLV,          // its logical volume
		  calorLV,          // its mother
		  kZAxis,           // axis of replication
		  fNofLayers,       // number of replica
		  layersThickness); // witdth of replica

  //                               
  // Lay1
  //
  G4VSolid* lay1S = new G4Box("Lay1",            // its name
			      calorSizeXY/2, calorSizeXY/2, lay1Thickness/2); // its size
                         
  G4LogicalVolume* lay1LV = new G4LogicalVolume(lay1S,        // its solid
						lay1Material, // its material
						"Lay1LV");        // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., -3.5*lay2Thickness), // its position
		    lay1LV,       // its logical volume                         
		    "Lay1",           // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  //                               
  // Lay2
  //
  G4VSolid* lay2S = new G4Box("Lay2",             // its name
			      calorSizeXY/2, calorSizeXY/2, lay2Thickness/2); // its size
                         
  G4LogicalVolume* lay2LV = new G4LogicalVolume(lay2S,             // its solid
						lay2Material,      // its material
						"Lay2LV");         // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., -2.5*lay1Thickness), // its position
		    lay2LV,            // its logical volume                         
		    "Lay2",            // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  //
  // Lay3
  //
  G4VSolid* lay3S = new G4Box("Lay3",             // its name
			      calorSizeXY/2, calorSizeXY/2, lay3Thickness/2); // its size
                         
  G4LogicalVolume* lay3LV = new G4LogicalVolume(lay3S,             // its solid
						lay3Material,      // its material
						"Lay3LV");         // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., -1.5*lay1Thickness), // its position
		    lay3LV,            // its logical volume                         
		    "Lay3",            // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  //
  // Lay4
  //
  G4VSolid* lay4S = new G4Box("Lay4",             // its name
			      calorSizeXY/2, calorSizeXY/2, lay4Thickness/2); // its size
                         
  G4LogicalVolume* lay4LV = new G4LogicalVolume(lay4S,             // its solid
						lay4Material,      // its material
						"Lay4LV");         // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., -.5*lay1Thickness), // its position
		    lay4LV,            // its logical volume                         
		    "Lay4",            // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  //
  // Lay5
  //
  G4VSolid* lay5S = new G4Box("Lay5",             // its name
			      calorSizeXY/2, calorSizeXY/2, lay5Thickness/2); // its size
                         
  G4LogicalVolume* lay5LV = new G4LogicalVolume(lay5S,             // its solid
						lay5Material,      // its material
						"Lay5LV");         // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., .5*lay1Thickness), // its position
		    lay5LV,            // its logical volume                         
		    "Lay5",            // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  //
  // Lay6
  //
  G4VSolid* lay6S = new G4Box("Lay6",             // its name
			      calorSizeXY/2, calorSizeXY/2, lay6Thickness/2); // its size
                         
  G4LogicalVolume* lay6LV = new G4LogicalVolume(lay6S,             // its solid
						lay6Material,      // its material
						"Lay6LV");         // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., 1.5*lay1Thickness), // its position
		    lay6LV,            // its logical volume                         
		    "Lay6",            // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  //
  // Lay7
  //
  G4VSolid* lay7S = new G4Box("Lay7",             // its name
			      calorSizeXY/2, calorSizeXY/2, lay7Thickness/2); // its size
                         
  G4LogicalVolume* lay7LV = new G4LogicalVolume(lay7S,             // its solid
						lay7Material,      // its material
						"Lay7LV");         // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., 2.5*lay1Thickness), // its position
		    lay7LV,            // its logical volume                         
		    "Lay7",            // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  //
  // Lay7
  //
  G4VSolid* lay8S = new G4Box("Lay8",             // its name
			      calorSizeXY/2, calorSizeXY/2, lay8Thickness/2); // its size
                         
  G4LogicalVolume* lay8LV = new G4LogicalVolume(lay8S,             // its solid
						lay8Material,      // its material
						"Lay8LV");         // its name
                                   
  new G4PVPlacement(0,                // no rotation
		    G4ThreeVector(0., 0., 3.5*lay1Thickness), // its position
		    lay8LV,            // its logical volume                         
		    "Lay8",            // its name
		    layerLV,          // its mother  volume
		    false,            // no boolean operation
		    0,                // copy number
		    fCheckOverlaps);  // checking overlaps 

  
  //
  // print parameters
  //
  /*
    G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << lay1Thickness/mm << "mm of " << lay1Material->GetName() 
    << " + "
    << lay2Thickness/mm << "mm of " << lay2Material->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  */
  //                                        
  // Visualization attributes
  //
  //  worldLV->SetVisAttributes(G4VisAttributes::Invisible);
  worldLV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,1.0)));

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  //  calorLV->SetVisAttributes(simpleBoxVisAtt);
  lay1LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  lay2LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  lay3LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  lay4LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  lay5LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  lay6LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  lay7LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  lay8LV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
  
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  CalorimeterSD* lay1SD = new CalorimeterSD("Lay1SD", "Lay1HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay1LV",lay1SD);

  CalorimeterSD* lay2SD = new CalorimeterSD("Lay2SD", "Lay2HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay2LV",lay2SD);

  CalorimeterSD* lay3SD = new CalorimeterSD("Lay3SD", "Lay3HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay3LV",lay3SD);

  CalorimeterSD* lay4SD = new CalorimeterSD("Lay4SD", "Lay4HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay4LV",lay4SD);

  CalorimeterSD* lay5SD = new CalorimeterSD("Lay5SD", "Lay5HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay5LV",lay5SD);

  CalorimeterSD* lay6SD = new CalorimeterSD("Lay6SD", "Lay6HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay6LV",lay6SD);

  CalorimeterSD* lay7SD = new CalorimeterSD("Lay7SD", "Lay7HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay7LV",lay7SD);

  CalorimeterSD* lay8SD = new CalorimeterSD("Lay8SD", "Lay8HitsCollection", fNofLayers);
  SetSensitiveDetector("Lay8LV",lay8SD);

  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
