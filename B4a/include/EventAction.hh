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
/// \file B4/B4a/include/EventAction.hh
/// \brief Definition of the B4a::EventAction class
#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4SDParticleFilter.hh" //a�adido por mi

namespace B4a
{

    /// Event action class
    ///
    /// It defines data members to hold the energy deposit and track lengths
    /// of charged particles in Absober and Gap layers:

    class EventAction : public G4UserEventAction
    {
    public:
        EventAction() = default;
        ~EventAction() override = default;

        void  BeginOfEventAction(const G4Event* event) override;
        void    EndOfEventAction(const G4Event* event) override;

        void AddDetector(G4double de, G4double dl);
        void SetEffectiveRange(G4double zmax);

    private:
        G4double  fEnergyDetector = 0.;
        G4double  fTrackLDetector = 0.;
        G4SDParticleFilter* neutronFilter;
        G4double fEffectiveRange = 0.;
    };

    // inline functions
    inline void EventAction::AddDetector(G4double de, G4double dl) {
        fEnergyDetector += de;
        fTrackLDetector += dl;
    }

    inline void EventAction::SetEffectiveRange(G4double zmax) {
        fEffectiveRange = zmax;
    }


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

