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
//
/// \file B4/B4a/src/DetectorConstruction.cc
/// \brief Implementation of the B4::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh" //Este lo he puesto yo
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh" //Este lo he añadido yo
#include "G4Transform3D.hh" //Este lo he añadido yo
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"


namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

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
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_B"); 
  nistManager->FindOrBuildMaterial("G4_AIR"); 
  nistManager->FindOrBuildMaterial("G4_Ge"); 
  nistManager->FindOrBuildMaterial("G4_SODIUM_IODIDE"); 
  nistManager->FindOrBuildMaterial("G4_CONCRETE"); 
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildElement("H");
  nistManager->FindOrBuildElement("C");
  nistManager->FindOrBuildElement("Cl");
  nistManager->FindOrBuildElement("O");
  nistManager->FindOrBuildElement("Mn"); //Manganeso
  nistManager->FindOrBuildElement("S"); //Azufre
  nistManager->FindOrBuildElement("Si"); //Silicio
  nistManager->FindOrBuildElement("Fe");
  nistManager->FindOrBuildElement("B");
  nistManager->FindOrBuildElement("Na");
  nistManager->FindOrBuildElement("I");
  nistManager->FindOrBuildElement("Al");
  nistManager->FindOrBuildElement("Ca");
  nistManager->FindOrBuildElement("K");
  nistManager->FindOrBuildElement("Mg");

 
  auto H = G4Element::GetElement("H");
  auto C = G4Element::GetElement("C");
  auto Cl = G4Element::GetElement("Cl");
  auto O = G4Element::GetElement("O");
  auto Mn = G4Element::GetElement("Mn");
  auto S = G4Element::GetElement("S");
  auto Si = G4Element::GetElement("Si");
  auto Fe = G4Element::GetElement("Fe");
  auto B = G4Element::GetElement("B");
  auto Na = G4Element::GetElement("Na");
  auto I = G4Element::GetElement("I");
  auto Mg = G4Element::GetElement("Mg");
  auto K = G4Element::GetElement("K");
  auto Ca = G4Element::GetElement("Ca");
  auto Al = G4Element::GetElement("Al");



  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  //new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);


  //PVC: C2H3Cl
  density = 1.4 * g / cm3;
  G4int ncomponents = 3;
  G4Material* PVC = new G4Material("PVC", density, ncomponents);
  PVC->AddElement(C, 2);
  PVC->AddElement(H, 3);
  PVC->AddElement(Cl, 1);

  //H2O
  
  density = 1.0 * g / cm3;
  G4int components = 2;
  G4Material* Water = new G4Material("Water", density, components); 
  Water->AddElement(H, 2);
  Water->AddElement(O, 1); 
  

  //Acero al carbono
  density = 7.85 * g / cm3;
  G4int numcomponents = 5;
  G4Material* steel = new G4Material("Steel", density, numcomponents);
  steel->AddElement(Fe, 98.0 * perCent);
  steel->AddElement(C, 0.2 * perCent);
  steel->AddElement(Mn, 1 * perCent);
  steel->AddElement(Si, 0.3 * perCent);
  steel->AddElement(S, 0.5 * perCent);


  //Acero al boro
  G4double density_b = 7.8 * g / cm3; // Densidad típica del acero
  G4Material* BoratedSteel = new G4Material("BoratedSteel", density_b, 3);
  BoratedSteel->AddElement(Fe, 98 * perCent);
  BoratedSteel->AddElement(C, 1 * perCent);
  BoratedSteel->AddElement(B, 1 * perCent);

  //NaI
      // Densidad del NaI
  G4double dens = 3.67 * g / cm3;

  // Crear el material NaI
  G4Material* NaI = new G4Material("NaI", dens, 2);
  NaI->AddElement(Na, 1);  // 1 átomo de sodio
  NaI->AddElement(I, 1);   // 1 átomo de iodo



  // Definir el material 'Cemento' (aproximado)
  G4double densit = 2.3 * g / cm3;  // Densidad típica del cemento
  G4Material* cement = new G4Material("Cement", densit, 9);

  // Composición típica del cemento Portland en fracciones de masa
  cement->AddElement(O, 0.521);  // 52.1% Oxígeno
  cement->AddElement(Si, 0.195);  // 19.5% Silicio
  cement->AddElement(Ca, 0.145);  // 14.5% Calcio
  cement->AddElement(Al, 0.043);  // 4.3% Aluminio
  cement->AddElement(Fe, 0.021);  // 2.1% Hierro
  cement->AddElement(Mg, 0.017);  // 1.7% Magnesio
  cement->AddElement(S, 0.015);  // 1.5% Azufre
  cement->AddElement(Na, 0.012);  // 1.2% Sodio
  cement->AddElement(K, 0.011);  // 1.1% Potasio



  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Ge");
  auto detectorMaterial = G4Material::GetMaterial("G4_SODIUM_IODIDE");
  auto air = G4Material::GetMaterial("G4_AIR");
  auto hormigon = G4Material::GetMaterial("G4_CONCRETE");
  auto iron = G4Material::GetMaterial("G4_Fe");
  auto PVC = G4Material::GetMaterial("PVC");
  auto Water = G4Material::GetMaterial("Water");
  auto Steel = G4Material::GetMaterial("Steel");
  auto BoratedSteel = G4Material::GetMaterial("BoratedSteel");
  auto NaI = G4Material::GetMaterial("NaI");
  auto Cement = G4Material::GetMaterial("Cement");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! detectorMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }
    

  //
  // World
  //
  G4double world_hx = 0.40 * m; 
  G4double world_hy = 0.40 * m;
  G4double world_hz = 0.40 * m;

  auto worldBox
      = new G4Box("World",           // its name
          world_hx, world_hy, world_hz); // its size

  auto worldLog
      = new G4LogicalVolume(
          worldBox,           // its solid
          defaultMaterial,  // its material
          "World");         // its name


  auto worldPhys = new G4PVPlacement(nullptr,  // no rotation
      G4ThreeVector(),                         // at (0,0,0)
      worldLog,                                 // its logical volume
      "World",                                 // its name
      nullptr,                                 // its mother  volume
      false,                                   // no boolean operation
      0,                                       // copy number
      fCheckOverlaps);                         // checking overlaps




  // Creo el detector:
  // 
  G4double detector_hx = 100.0 * cm;
  G4double detector_hy = 110.0 * cm;
  G4double detector_hz = 100.0 * cm;

  G4Box* detectorBox = new G4Box("Detector", detector_hx, detector_hy, detector_hz);
  G4LogicalVolume* detectorLog
      = new G4LogicalVolume(detectorBox, Cement, "Detector");


  //Coloco el detector justo detrás de el cilindro
  G4double detectorPos_x = 0.0 * m; 
  G4double detectorPos_y = 0.0 * m;
  G4double detectorPos_z = +50.0 * cm;   
    
  //Volumen fisico del detector:
  //
      fDetectorPhys
      = new G4PVPlacement(0,                       // no rotation
          G4ThreeVector(detectorPos_x, detectorPos_y, detectorPos_z),
          // translation position
          detectorLog,              // its logical volume
          "Detector",               // its name
          worldLog,                // its mother (logical) volume
          false,                   // no boolean operations
          0);                      // its copy number




    // Visualization attributes
    //G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red
    //worldVisAtt->SetVisibility(true);
    //worldLog->SetVisAttributes(worldVisAtt);
  worldLog->SetVisAttributes(G4VisAttributes::GetInvisible());


    G4VisAttributes* detectorVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // White
    detectorVisAtt->SetVisibility(true);
    detectorLog->SetVisAttributes(detectorVisAtt);
    

    



  //
  // Always return the physical World
  //
   return worldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



}

