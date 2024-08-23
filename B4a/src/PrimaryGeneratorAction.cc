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
/// \file B4/B4a/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B4::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh" //incluida por mi
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace B4
{

    PrimaryGeneratorAction::PrimaryGeneratorAction()
    {
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        // default particle kinematic
        //
        auto particleDefinition
            = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        fParticleGun->SetParticleDefinition(particleDefinition);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

        // Establece la energía de los neutrones
        fParticleGun->SetParticleEnergy(2.5 * MeV);

        // Establece la posición inicial de la pistola de partículas
        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -25 * cm));
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    PrimaryGeneratorAction::~PrimaryGeneratorAction()
    {
        delete fParticleGun;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
    {

        G4double worldZHalfLength = 0.;
        auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");


        // Check that the world volume has box shape
        G4Box* worldBox = nullptr;
        if (worldLV) {
            worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
        }


        if (worldBox) {
            worldZHalfLength = worldBox->GetZHalfLength();
        }
        else {
            G4ExceptionDescription msg;
            msg << "World volume of box shape not found." << G4endl;
            msg << "Perhaps you have changed geometry." << G4endl;
            msg << "The gun will be place in the center.";
            G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
                "MyCode0002", JustWarning, msg);
        }

        G4double cosMin = std::cos(12.0 * deg); // cos(20 grados)
        G4double cosMax = 1.0;                 // cos(0 grados)

        // Generar un número aleatorio entre cosMin y cosMax
        G4double randomNumber = cosMin + (cosMax - cosMin) * G4UniformRand();

        // Calcular el ángulo correspondiente usando arccos (en radianes)
        G4double theta = std::acos(randomNumber);

        // Generar un ángulo phi aleatorio entre 0 y 360 grados
        G4double phi = G4UniformRand() * 360.0 * deg;

        // Convertir coordenadas esféricas a cartesianas
        G4double ux = std::sin(theta) * std::cos(phi);
        G4double uy = std::sin(theta) * std::sin(phi);
        G4double uz = std::cos(theta);

        // Establecer la dirección del momento de la partícula
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

        // Generate the primary vertex
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}
