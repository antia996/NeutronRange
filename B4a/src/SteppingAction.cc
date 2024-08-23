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
/// \file B4/B4a/src/SteppingAction.cc
/// \brief Implementation of the B4a::SteppingAction class
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4AnalysisManager.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include <set>

using namespace B4;

namespace B4a
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    SteppingAction::SteppingAction(const DetectorConstruction* detConstruction,
        EventAction* eventAction)
        : fDetConstruction(detConstruction),
        fEventAction(eventAction),
        range(0.0) //inicializa zmax en 0
    {}

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void SteppingAction::UserSteppingAction(const G4Step* step)
    {
        //Este metodo es la parte central del codigo. Se llama cada vez que una partícula da un "paso" en la simulacion
        // Collect energy and track length step by step
        
          // get volume of the current step
        auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(); //obtiene el volumen del paso actual
        auto nextVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume(); //obtiene el volumen del siguiente paso

        // energy deposit
        auto edep = step->GetTotalEnergyDeposit(); //obtiene la energía depositada en este paso

        auto kineticEnergy = step->GetPreStepPoint()->GetKineticEnergy(); //obtiene la energía cinética de la partícula antes del paso

        auto particle = step->GetTrack()->GetDefinition();

        G4int trackID = step->GetTrack()->GetTrackID();

        G4ThreeVector position = step->GetPostStepPoint()->GetPosition();

        //G4cout << "Track ID: " << trackID << " - Particle: " << particle->GetPDGEncoding() << " - Position post: " << position.x() << ":" << position.y() << ":" << position.z() << G4endl;
        //G4cout << "Position pre: " << step->GetPreStepPoint()->GetPosition().x() << ":" << step->GetPreStepPoint()->GetPosition().y() << ":" << step->GetPreStepPoint()->GetPosition().z() << G4endl;
        // step length
        
        G4double stepLength = 0.;
        //if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0. && particle->GetPDGEncoding() == 2112) {
        if (trackID == 1 && volume == fDetConstruction->GetDetectorPhys()) {
            // En el primer paso de la partícula primaria, reiniciar `range`
            if (step->GetTrack()->GetCurrentStepNumber() == 1) {
                range = position.z();
            }


            stepLength = step->GetStepLength();

            fEventAction->AddDetector(edep, stepLength);

            //G4cout << position.x() << ":" << position.y() << ":" << position.z() << G4endl;

            if (position.z() > range) range = position.z();
            // Only consider neutrons

                    // Verifica si la partícula ha finalizado su trayecto
            if (step->GetTrack()->GetTrackStatus() == fStopAndKill ||
                step->GetTrack()->GetTrackStatus() == fStopButAlive ||
                step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != fDetConstruction->GetDetectorPhys()){
                //nextVolume != fDetConstruction->GetDetectorPhys()) {
                
                // Imprime la posición máxima en Z al final de la trayectoria
                G4cout << "Max Z Position for trackID " << trackID << ": " << range << G4endl;
                fEventAction->SetEffectiveRange(range);
            }

            
            /*
            if (volume == fDetConstruction->GetDetectorPhys()) {
                fEventAction->AddDetector(edep, stepLength);
                //G4cout << "Antes: " << fPreviousPositionMap[trackID].z() << G4endl;
                if (fInitialPositionMap.find(trackID) == fInitialPositionMap.end()) { //se almacena la posicion inicial del neutron si no esta registrada
                    // If this track ID is not in the map, store the initial position
                    fInitialPositionMap[trackID] = step->GetPreStepPoint()->GetPosition();
                }

                // Track position to fill H2 histogram
                G4double x = step->GetPreStepPoint()->GetPosition().x();
                G4double y = step->GetPreStepPoint()->GetPosition().y();
                auto analysisManager = G4AnalysisManager::Instance();
                analysisManager->FillH2(0, x, y);

                // Almacena la posición actual como la posición anterior para el siguiente paso
                fPreviousPositionMap[trackID] = step->GetPreStepPoint()->GetPosition();
                //G4cout << "Despues: " << fPreviousPositionMap[trackID].z() << G4endl;

            }
            */
            

            //if (fInitialPositionMap.find(trackID) != fInitialPositionMap.end()) {
                // If the particle exits the detector, calculate the effective range
                //G4double initialPosition = 0.0;
                //G4double initialPosition = fInitialPositionMap[trackID].z();
                //G4double previousPosition = fPreviousPositionMap[trackID].z();
                //G4double finalPosition = step->GetPreStepPoint()->GetPosition().z();
                //G4double range = std::abs(finalPosition - initialPosition);
                //G4cout << finalPosition << "\t";
                //fEventAction->AddEffectiveRange(range);
                //G4cout << "initialPosition " << initialPosition << ", finalPosition: " << previousPosition << ", range: " << range << G4endl;
                

                //if (finalPosition > range) {
                    //range = finalPosition;
                    //G4double range = zmax;
                    //G4cout << "la z maxima es:" << zmax << G4endl;
                    
                //}

                //if (fPrintedTrackIDs.find(trackID) == fPrintedTrackIDs.end()) {
                    //G4cout << "La z maxima es: " << range << G4endl;
                    //fPrintedTrackIDs.insert(trackID);
                    //fEventAction->SetEffectiveRange(range);
                //}
                // 
                //if (fInitialPositionMap.find(trackID) != fInitialPositionMap.end()) {
                    //G4cout << "La z maxima es: " << range << G4endl;
                    //fEventAction->SetEffectiveRange(range);
                //}

                /*
                if (zmax > range) {
                    range = zmax;
                    //G4cout << "zmax:" << range << G4endl;
                    if (fPrintedTrackIDs.find(trackID) == fPrintedTrackIDs.end()) {
                        G4cout << "La z maxima es: " << range << G4endl;
                        fPrintedTrackIDs.insert(trackID);
                        fEventAction->SetEffectiveRange(range);
                    }
                }
                */
                
              
                
                /*
                
                //Compruebo si el track ha finalizado , es decir, si la partícula ha salido del detector o se ha detenido
                if (!nextVolume || nextVolume != fDetConstruction->GetDetectorPhys()) {
                    G4cout << "La mayor Z es: " << zmax << G4endl;
                    zmax = 0.0;
                    fInitialPositionMap.erase(trackID);
                    fPreviousPositionMap.erase(trackID);
                }
                */
                /*
                // Bucle para encontrar el rango máximo mientras el rango sea mayor que 0
                while (range > 0) {
                    if (range > maxRange) {
                        maxRange = range;
                        G4cout << "el rango maximo es: " << maxRange << G4endl;
                        //fEventAction->AddEffectiveRange(maxRange);
                    }
                    break;  
                }
                */
                // Elimina el ID de la trayectoria del mapa

                //fInitialPositionMap.erase(trackID);
                //fPreviousPositionMap.erase(trackID);
            //}
            
        }
    
        //fEventAction->SetEffectiveRange(range);
        //G4cout << "La mayor Z es:" << range << G4endl;
        
        
    } 
    

    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


