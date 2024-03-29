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
/// \file PicoVeto/src/PicoVetoStackingAction.cc
/// \brief Implementation of the PicoVetoStackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PicoVetoStackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "PicoVetoSensitiveDetector.hh"

G4int PicoVetoStackingAction::g_N_optical_trajectories = 0; //zero of global # of optical photons

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoStackingAction::PicoVetoStackingAction()
  : G4UserStackingAction(),
    fScintillationCounter(0), fCerenkovCounter(0),  fTotalCerenkovCounter(0), fCutTracksCounter(0)

{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoStackingAction::~PicoVetoStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
PicoVetoStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  //  std::ofstream log("logfile.txt", std::ios_base::app | std::ios_base::out);
  //  if(aTrack->GetDefinition() == G4MuonMinus::MuonMinusDefinition())
  //  { log << aTrack->GetPosition().x() << "\t" <<aTrack->GetPosition().y() << "\t" << aTrack->GetPosition().z()<< "\n"; }


  fCutTracksCounter++;



  
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
    { // particle is optical photon
      if(aTrack->GetParentID()>0)
	{ // particle is secondary
	  if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
	    fScintillationCounter++;
	  if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
	    {


	           if( fCutTracksCounter%1000 != 0)
	          	return fKill;
	          else{

		fCerenkovCounter++;
		fTotalCerenkovCounter++;
		g_N_optical_trajectories++;}

	    } 
	}
    }



    
  
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PicoVetoStackingAction::NewStage()
{
  //G4cout << "Number of Scintillation photons produced in this event : "
  //<< fScintillationCounter << G4endl;
  G4cout << G4endl;
  G4cout << "Number of Cerenkov photons produced in this event : "
         << fCerenkovCounter << G4endl;
 G4cout << "total number of Cerenkov photons produced until now : "
         << fTotalCerenkovCounter << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PicoVetoStackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
