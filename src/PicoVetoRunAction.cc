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
/// \file PicoVeto/src/PicoVetoRunAction.cc
/// \brief Implementation of the PicoVetoRunAction class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make this appear first!
#include "G4Timer.hh"

#include "PicoVetoRunAction.hh"

#include "G4Run.hh"
#include "G4ScoringManager.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include <assert.h>
#include "PicoVetoSensitiveDetector.hh"
#include <time.h>   
G4int PicoVetoRunAction::g_run_N = 0; //zero of run number
//G4Timer* PicoVetoRunAction::myTimer ; 
time_t PicoVetoRunAction::myTimer ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoRunAction::PicoVetoRunAction()
 : G4UserRunAction(),
   fTimer(0)
{
  fTimer = new G4Timer;
 // myTimer = new G4Timer;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoRunAction::~PicoVetoRunAction()
{
  delete fTimer, myTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void PicoVetoRunAction::BeginOfRunAction(const G4Run* aRun)
{
  fTimer->Start();  
 // myTimer->Start();

  //  G4cout << "fTimer : " << *fTimer <<  G4endl;

    time(&myTimer);  /* get current time; same as: timer = time(NULL)  */

   // G4cout << " Clock UTC : " << myTimer  << G4endl;
	//G4cout << " Clock starts at : " << ctime(&myTimer)  << G4endl;

	  G4cout << "### Run " << aRun->GetRunID() << " start on  " <<  ctime(&myTimer)  << G4endl;

	
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PicoVetoRunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
 // myTimer->Stop();

  G4cout << "### Run  " << aRun->GetRunID() << " Finished! " << G4endl;
  G4cout << " Timer : " << *fTimer << G4endl << G4endl;
//	  G4cout << " Timer : " << *myTimer << G4endl << G4endl;
	 

  
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
