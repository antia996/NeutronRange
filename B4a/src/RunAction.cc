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
/// \file B4/B4a/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();
  

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  analysisManager->CreateH1("EdepDetector", "Edep in detector", 3000, 0., 15 * MeV);
  analysisManager->CreateH1("EdepTracker", "Edep in tracker", 3000, 0., 15 * MeV);
  analysisManager->CreateH1("LengthDetector", "trackL in detector", 3000, 0., 50 * cm);
  analysisManager->CreateH1("LengthTracker", "trackL in tracker", 3000, 0., 50 * cm);
  //analysisManager->CreateH1("EdepBackDetector", "Edep in back detector", 3000, 0., 10 * MeV);
  //analysisManager->CreateH1("LengthBackDetector", "trackL in back detector", 3000, 0., 50*cm);

  //Queremos crear un histograma de las posiciones X e Y de las particulas que llegan al detector
  analysisManager->CreateH2("h2", "Posiciones de las particulas en el detector", 1000., -150., 150., 1000., -150., 150.);
  //analysisManager->CreateH2("h2Back", "Posiciones de las partículas en el back detector", 1000., -150., 150., 1000, -150., 150.);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("B4", "Edep, TrackL and positions");
  analysisManager->CreateNtupleDColumn("EDetector");
  analysisManager->CreateNtupleDColumn("ETracker");
  analysisManager->CreateNtupleDColumn("LDetector");
  analysisManager->CreateNtupleDColumn("LTracker");
  //analysisManager->CreateNtupleDColumn("EBackDetector");
  //analysisManager->CreateNtupleDColumn("LBackDetector");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "B4.root";
  // Other supported output types:
  // G4String fileName = "B4.csv";
  // G4String fileName = "B4.hdf5";
  // G4String fileName = "B4.xml";
  analysisManager->OpenFile(fileName);
  G4cout << "Using " << analysisManager->GetType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1)) { //aqui he añadido analysisManager->GetH2(1) no se si esta bien
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl;
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl;
    }

    //Esto lo añado yo:
    G4cout << " EDetector : mean = "
        << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
        << " rms = "
        << G4BestUnit(analysisManager->GetH1(0)->rms(), "Energy") << G4endl;

    G4cout << " ETracker : mean = "
        << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
        << " rms = "
        << G4BestUnit(analysisManager->GetH1(1)->rms(), "Energy") << G4endl;

    G4cout << " LDetector : mean = "
        << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length")
        << " rms = "
        << G4BestUnit(analysisManager->GetH1(2)->rms(), "Length") << G4endl;

    G4cout << " LTracker : mean = "
        << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
        << " rms = "
        << G4BestUnit(analysisManager->GetH1(3)->rms(), "Length") << G4endl;    
    //G4cout << " EBackDetector : mean = "
      //  << G4BestUnit(analysisManager->GetH1(4)->mean(), "Energy")
       // << " rms = "
        //<< G4BestUnit(analysisManager->GetH1(4)->rms(), "Energy") << G4endl;
    //G4cout << " LBackDetector : mean = "
      //  << G4BestUnit(analysisManager->GetH1(5)->mean(), "Length")
        //<< " rms = "
        //<< G4BestUnit(analysisManager->GetH1(5)->rms(), "Length") << G4endl;
  }
  

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
