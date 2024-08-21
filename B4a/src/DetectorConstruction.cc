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
 
  auto H = G4Element::GetElement("H");
  auto C = G4Element::GetElement("C");
  auto Cl = G4Element::GetElement("Cl");
  auto O = G4Element::GetElement("O");
  auto Mn = G4Element::GetElement("Mn");
  auto S = G4Element::GetElement("S");
  auto Si = G4Element::GetElement("Si");
  auto Fe = G4Element::GetElement("Fe");
  auto B = G4Element::GetElement("B");


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

  if ( ! defaultMaterial || ! absorberMaterial || ! detectorMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }
    

  //
  // World
  //
  G4double world_hx = 0.50 * m; 
  G4double world_hy = 0.50 * m;
  G4double world_hz = 0.50 * m;

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


  //Creo la columna:
  // 
  G4double tracker_hx = 5 * cm;
  G4double tracker_hy = 10 * cm;
  G4double tracker_hz = 5 * cm;

  auto trackerTube
      = new G4Box("Tracker", 
          tracker_hx, tracker_hy, tracker_hz); 

  auto trackerTubeLV = new G4LogicalVolume(trackerTube, hormigon, "trackerTubeLV");

  auto trackerTubePV
      = new G4PVPlacement(0,
          G4ThreeVector(),
          trackerTubeLV,
          "trackerTubePV",
          worldLog,
          false,
          0,
          fCheckOverlaps);

  //Voy a recortar esta columna de hormigon 
  G4double initRadius = 0.0 * m;
  G4double finRadius = 6.0 * mm;
  G4double high = tracker_hy;
  G4double initAngle = 0. * deg;
  G4double finAngle = 360. * deg;
  
  //Solido:
  auto SteelTube
      = new G4Tubs("SteelTube",
          initRadius,
          finRadius,
          high,
          initAngle,
          finAngle);

  auto SteelTubeLV = new G4LogicalVolume(SteelTube, BoratedSteel, "SteelTubeLV");

  G4RotationMatrix* rot90 = new G4RotationMatrix();
  rot90->rotateX(90.0 * deg);  // Girar 90 grados alrededor del eje Y  
  auto SteelTubePV
      = new G4PVPlacement(rot90,
          G4ThreeVector(),
          SteelTubeLV,
          "SteelTubeLV",
          trackerTubeLV,
          false,
          0,
          fCheckOverlaps);


  

  //Hago otro agujero en la columna
  G4double R_0 = 0.0 * m;
  G4double R_f = 6 * mm;
  G4double high_0 = tracker_hy;
  G4double Angle_0 = 0. * deg;
  G4double Angle_fin = 360 * deg;

  auto AirTube
      = new G4Tubs("AirTube",
          R_0,
          R_f,
          high_0,
          Angle_0,
          Angle_fin);

  auto AirTubeLV = new G4LogicalVolume(AirTube, air, "AirTubeLV");

  auto AirTubePV = new G4PVPlacement(rot90,
      G4ThreeVector(-3.5*cm, 0.0, 0.0),
      AirTubeLV,
      "AirTubePV",
      trackerTubeLV,
      false,
      0,
      fCheckOverlaps);


  G4double Rinit = 0.0;
  G4double Rfin = 6 * mm;
  G4double altura = tracker_hy;
  G4double Angleinit = 0. * deg;
  G4double Anglefin = 360. * deg;

  auto WaterTube
      = new G4Tubs("WaterTube",
          Rinit,
          Rfin,
          altura,
          Angleinit,
          Anglefin);

  auto WaterTubeLV = new G4LogicalVolume(WaterTube, Water, "WaterTubeLV");

  auto WaterTubePV
      = new G4PVPlacement(rot90,
          G4ThreeVector(+3.5*cm, 0.0, 0.0),
          WaterTubeLV,
          "WaterTubeLV",
          trackerTubeLV,
          false,
          0,
          fCheckOverlaps);



          
  /*

  //Voy a recortar esta columna de hormigon 
  G4double initRadius = 0.0 * m;
  G4double finRadius = 12 * mm;
  G4double high = tracker_hy;
  G4double initAngle = 0. * deg;
  G4double finAngle = 360. * deg;


  //Solido:
  auto AirTube
      = new G4Tubs("AirTube",
          initRadius,
          finRadius,
          high,
          initAngle,
          finAngle);



  //Hago otro agujero en la columna
  G4double R_0 = 0.0 * m; 
  G4double R_f = 12 * mm; 
  G4double high_0 = tracker_hy; 
  G4double Angle_0 = 0. * deg; 
  G4double Angle_fin = 360 * deg; 

  auto AirTube2
      = new G4Tubs("AirTube2",
          R_0,
          R_f,
          high_0,
          Angle_0,
          Angle_fin);

  //Hago un tercer agujero en la columna
  G4double Rinit = 0.0;
  G4double Rfin = 5.0 * cm;
  G4double h = tracker_hy;
  G4double Ang_0 = 0 * deg;
  G4double Ang_fin = 360 * deg;

  auto AirTube3
      = new G4Tubs("AirTube3",
          Rinit,
          Rfin,
          high_0,
          Angle_0,
          Angle_fin);




  G4RotationMatrix* rot90 = new G4RotationMatrix();
  rot90->rotateX(90.0 * deg);  // Girar 90 grados alrededor del eje Y  
  G4SubtractionSolid* trackerTubeWithHole = new G4SubtractionSolid("trackerTubeWithHole",
      trackerTube,
      AirTube,
      rot90,
      G4ThreeVector());

  G4double trackerPosX = 0.0 * cm;
  G4double trackerPosY = 0.0 * cm;
  G4double trackerPosZ = 0.0 * cm;

  G4SubtractionSolid* trackerTubeWith2Holes = new G4SubtractionSolid("trackerTubeWith2Holes",
      trackerTubeWithHole,
      AirTube2,
      rot90,
      G4ThreeVector(-15.0*cm, 0.0, 00.0*cm));

  G4SubtractionSolid* trackerTubeWith3Holes = new G4SubtractionSolid("trackerTubeWith3Holes",
      trackerTubeWith2Holes,
      AirTube3,
      rot90,
      G4ThreeVector(+15.0*cm, 0.0, 0.0));


  auto trackerLV = new G4LogicalVolume(trackerTubeWith3Holes, air, "trackerLV");

  G4PVPlacement* trackerPV = new G4PVPlacement(0,
      G4ThreeVector(),
      trackerLV,
      "trackerPV",
      worldLog,
      false,
      0,
      fCheckOverlaps);




  //Ahora voy a rellenar los huecos con barras de acero, aire y agua

  auto AirTube4
      = new G4Tubs("AirTube4",
          initRadius,
          finRadius,
          high,
          initAngle,
          finAngle);

  auto AirTube4LV = new G4LogicalVolume(AirTube4, hormigon, "AirTube4LV");

  auto AirTube4PV
      = new G4PVPlacement(rot90,
          G4ThreeVector(),
          AirTube4LV,
          "AirTube4PV",
          trackerLV,
          false,
          0,
          fCheckOverlaps);
  
  auto WaterTube
      = new G4Tubs("WaterTube",
          R_0,
          R_f,
          high_0,
          Angle_0,
          Angle_fin);

  auto WaterTubeLV = new G4LogicalVolume(WaterTube, Water, "WaterTubeLV");

  auto WaterTubePV = new G4PVPlacement(rot90,
      G4ThreeVector(-15*cm, 0.0,0.0),
      WaterTubeLV,
      "WaterTubePV",
      trackerLV,
      false,
      0,
      fCheckOverlaps);

  auto SteelTube
      = new G4Tubs("SteelTube",
          Rinit,
          Rfin,
          high_0,
          Angle_0,
          Angle_fin);

  auto SteelTubeLV = new G4LogicalVolume(SteelTube, Steel, "SteelTubeLV");

  auto SteelTubePV = new G4PVPlacement(rot90,
      G4ThreeVector(+15.0*cm, 0.0, 0.0),
      SteelTubeLV,
      "SteelTubeLV",
      trackerLV,
      false,
      0,
      fCheckOverlaps);




    */


 /*
  //Ahora creamos el volumen logico:
  // 
  auto trackerLog
      = new G4LogicalVolume(trackerTube, hormigon, "Tracker"); 


  // Creamos el volumen fisico
  // 
  G4double pos_x = 0.0*m;
  G4double pos_y = 0.0*m;
  G4double pos_z = 0.0*m;

  //G4RotationMatrix* rot90 = new G4RotationMatrix(); 
  //rot90->rotateX(90.0 * deg);  // Girar 90 grados alrededor del eje Y       
  
  ftrackerPhys = new G4PVPlacement(nullptr,  // no rotation
      G4ThreeVector(),          // at (0,0,0)
      trackerLog,                  // its logical volume
      "Tracker",            // its name
      worldLog,                  // its mother  volume
      false,                    // no boolean operation
      0,                        // copy number
      fCheckOverlaps);          // checking overlaps

  
  //Voy a crear una tubería de PVC o cobre, de 10 cm de diametro y 3 mm de perfil que pueda tener dentro agua o aire
  G4double initRadius = 0.0497 * m;
  G4double finRadius = 0.050 * m;
  G4double high = 0.5 * m;
  G4double initAngle = 0. * deg;
  G4double finAngle = 360. * deg;


  //Solido:
  auto trackerIron
      = new G4Tubs("Iron",
          initRadius,
          finRadius,
          high,
          initAngle,
          finAngle);

  //Ahora creamos el volumen logico:
  // 
  auto trackerIronLog
      = new G4LogicalVolume(trackerIron, PVC, "Iron"); 


  // Creamos el volumen fisico
  // 
  G4double xinit = 0.0 * m;
  G4double yinit = 0.0 * m;
  G4double zinit = 0.0 * m;
  G4RotationMatrix* rot90 = new G4RotationMatrix();
  rot90->rotateX(90.0 * deg);  // Girar 90 grados alrededor del eje Y        

  ftrackerIronPhys = new G4PVPlacement(rot90,  // no rotation 
      G4ThreeVector(),          // at (0,0,0)
      trackerIronLog,                  // its logical volume
      "Iron",            // its name
      worldLog,                  // its mother  volume
      false,                    // no boolean operation
      0,                        // copy number
      fCheckOverlaps);          // checking overlaps






  //dentro de este volumen quiero que haya aire o agua, así que hago un cilindro macizo de aire o agua

  G4double init_Radius = 0. * cm;
  G4double fin_Radius = 4.96 * cm;
  G4double high_h = 0.5 * m; 
  G4double init_Angle = 0. * deg;
  G4double fin_Angle = 360. * deg;


  //Solido:
  auto trackerWaterOrAir
      = new G4Tubs("Water or Air",
          init_Radius,
          fin_Radius,
          high_h, 
          init_Angle,
          fin_Angle);

  //Ahora creamos el volumen logico:
  // 
  auto trackerWaterOrAirLog
      = new G4LogicalVolume(trackerWaterOrAir, Water, "Water or Air");  


  // Creamos el volumen fisico
  // 
  G4double x_init = 0.0 * m;
  G4double y_init = 0.0 * m;
  G4double z_init = 0.0 * m;
  G4RotationMatrix* matrixrot90 = new G4RotationMatrix();
  matrixrot90->rotateX(90.0 * deg);  // Girar 90 grados alrededor del eje Y        

  auto ftrackerWaterOrAirPhys = new G4PVPlacement(matrixrot90,  // no rotation 
      G4ThreeVector(),          // at (0,0,0)
      trackerWaterOrAirLog,                  // its logical volume
      "Water or Air",            // its name
      trackerLog,                  // its mother  volume
      false,                    // no boolean operation
      0,                        // copy number
      fCheckOverlaps);          // checking overlaps


*/





  // Creo el detector:
  // 
  G4double detector_hx = 7.0 * cm;
  G4double detector_hy = 15 * cm;
  G4double detector_hz = 7.0 * cm;

  G4Box* detectorBox = new G4Box("Detector", detector_hx, detector_hy, detector_hz);
  G4LogicalVolume* detectorLog
      = new G4LogicalVolume(detectorBox, detectorMaterial, "Detector");


  //Coloco el detector justo detrás de el cilindro
  G4double detectorPos_x = 0.0 * m; 
  G4double detectorPos_y = 0.0 * m;
  G4double detectorPos_z = tracker_hz + 10.0*cm;   
    
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

/*

// Creo OTRO DETECTOR: BACK DETECTOR:
  // 
  G4double BackDetector_hx = 7.0 * cm;
  G4double BackDetector_hy = 15 * cm;
  G4double BackDetector_hz = 7.0 * cm;

  G4Box* BackDetectorBox = new G4Box("BackDetector", BackDetector_hx, BackDetector_hy, BackDetector_hz);
  G4LogicalVolume* BackDetectorLog 
      = new G4LogicalVolume(detectorBox, detectorMaterial, "BackDetector");


  //Coloco el detector justo detrás de el cilindro
  G4double BackDetectorPos_x = 0.0 * m; 
  G4double BackDetectorPos_y = 0.0 * m;
  G4double BackDetectorPos_z = tracker_hz + 30 * cm;   //con 0.10 m están pegados 
    
  //Volumen fisico del detector:
  //
      fBackDetectorPhys
      = new G4PVPlacement(0,                       // no rotation
          G4ThreeVector(BackDetectorPos_x, BackDetectorPos_y, BackDetectorPos_z),
          // translation position
          BackDetectorLog,              // its logical volume
          "BackDetector",               // its name
          worldLog,                // its mother (logical) volume
          false,                   // no boolean operations
          0);                      // its copy number




//Voy a introducir un pequeño cilindro de hierro dentro de la columna de hormigon:

  G4double initRadius = 0. * m; 
  G4double finRadius = 2.5 * cm;
  G4double high = 0.5 * m; 
  G4double initAngle = 0. * deg; 
  G4double finAngle = 360. * deg;


//Solido:
  auto trackerIron
      = new G4Tubs("Iron",
          initRadius,
          finRadius, 
          high,
          initAngle, 
          finAngle); 

//Ahora creamos el volumen logico:
// 
  auto trackerIronLog
      = new G4LogicalVolume(trackerIron, iron, "Iron");


// Creamos el volumen fisico
// 
  G4double xinit = 0.0 * m; 
  G4double yinit = 0.0 * m; 
  G4double zinit = 0.0 * m; 

  G4RotationMatrix* matrixrot90 = new G4RotationMatrix(); 
  matrixrot90->rotateX(90.0 * deg);  // Girar 90 grados alrededor del eje Y        

  ftrackerIronPhys = new G4PVPlacement(matrixrot90,  // no rotation 
      G4ThreeVector(),          // at (0,0,0)
      trackerIronLog,                  // its logical volume
      "Iron",            // its name
      worldLog,                  // its mother  volume
      false,                    // no boolean operation
      0,                        // copy number
      fCheckOverlaps);          // checking overlaps
*/




  //
  // print parameters
  //
  //G4cout
    //<< G4endl
    //<< "------------------------------------------------------------" << G4endl
    //<< "---> The calorimeter is " << nofLayers << " layers of: [ "
    //<< absoThickness/mm << "mm of " << absorberMaterial->GetName()
    //<< " + "
    //<< gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    //<< "------------------------------------------------------------" << G4endl;


    // Visualization attributes
    //G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red
    //worldVisAtt->SetVisibility(true);
    //worldLog->SetVisAttributes(worldVisAtt);
  worldLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  /*
    G4VisAttributes* trackerVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));// Green
    trackerVisAtt->SetVisibility(true);
    trackerLog->SetVisAttributes(trackerVisAtt);
*/

    G4VisAttributes* detectorVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // White
    detectorVisAtt->SetVisibility(true);
    detectorLog->SetVisAttributes(detectorVisAtt);
    
    /*
    G4VisAttributes* cilVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); //Azul
    cilVisAtt->SetVisibility(true);
    trackerWaterOrAirLog->SetVisAttributes(cilVisAtt);
    

    G4VisAttributes* BackDetectorVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    BackDetectorVisAtt->SetVisibility(true);
    BackDetectorLog->SetVisAttributes(BackDetectorVisAtt);
    */
  
    G4VisAttributes* WaterTubeVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
    WaterTubeVisAtt->SetVisibility(true);
    WaterTubeLV->SetVisAttributes(WaterTubeVisAtt);

    G4VisAttributes* SteelTubeVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    SteelTubeVisAtt->SetVisibility(true);
    SteelTubeLV->SetVisAttributes(SteelTubeVisAtt);

    G4VisAttributes* AirTube4VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    AirTube4VisAtt->SetVisibility(true);
    AirTubeLV->SetVisAttributes(AirTube4VisAtt);

    G4VisAttributes* trackerVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    trackerVisAtt->SetVisibility(true);
    trackerTubeLV->SetVisAttributes(trackerVisAtt);
 



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

