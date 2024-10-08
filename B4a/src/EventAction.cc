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
/// \file B4/B4a/src/EventAction.cc
/// \brief Implementation of the B4a::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <map>

namespace B4a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
  // initialisation per event
  fEnergyDetector = 0.;
  fTrackLDetector = 0.;
  fEffectiveRange = 0.;
  G4cout << "Evento: " << event->GetEventID() << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{


    G4int nPrimaries = event->GetNumberOfPrimaryVertex();

    for (G4int iVertex = 0; iVertex < nPrimaries; ++iVertex) {
        G4PrimaryVertex* vertex = event->GetPrimaryVertex(iVertex);

        G4PrimaryParticle* primary = vertex->GetPrimary();

        
        //G4cout << "PDG code: " << primary->GetPDGcode() << G4endl;

        if (primary->GetPDGcode() == 2112){
            // get analysis manager
                auto analysisManager = G4AnalysisManager::Instance();
                
            // fill histograms
            analysisManager->FillH1(0, fEnergyDetector);
            analysisManager->FillH1(1, fTrackLDetector);
            analysisManager->FillH1(2, fEffectiveRange);

            // fill ntuple
            analysisManager->FillNtupleDColumn(0, fEnergyDetector);
            analysisManager->FillNtupleDColumn(1, fTrackLDetector);
            analysisManager->FillNtupleDColumn(2, fEffectiveRange);
            analysisManager->AddNtupleRow();

            // Print per event (modulo n)
            //
            auto eventID = event->GetEventID();
            auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();

            
            /*if ((printModulo > 0) && (eventID % printModulo == 0)) {
                G4cout
                    << "   Detector: total energy: " << std::setw(7)
                    << G4BestUnit(fEnergyDetector, "Energy")
                    << "       total track length: " << std::setw(7)
                    << G4BestUnit(fTrackLDetector, "Length")
                    << G4endl              
                    << G4endl;

                G4cout << "--> End of event " << eventID << "\n" << G4endl;
            }*/


        }
        // Debug output
        //G4cout << "Effective Range: " << fEffectiveRange << G4endl;
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

} 
