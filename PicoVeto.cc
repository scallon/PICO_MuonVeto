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
/// \file PicoVeto/PicoVeto.cc
/// \brief Main program of the PicoVeto example
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Description: Test of Continuous Process G4Cerenkov
//              and RestDiscrete Process G4Scintillation
//              -- Generation Cerenkov Photons --
//              -- Generation Scintillation Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     1996-04-30
// Author:      Juliet Armstrong
// mail:        gum@triumf.ca
//     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "PicoVetoPhysicsList.hh"
#include "PicoVetoDetectorConstruction.hh"

#include "PicoVetoActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4ScoringManager.hh"


#include "PicoVetoSensitiveDetector.hh"
#include "PicoVetoEventAction.hh"

#include "PicoVetoPrimaryGeneratorAction.hh"
#include "PicoVetoRunAction.hh"
#include "PicoVetoStackingAction.hh"
#include "PicoVetoSteppingVerbose.hh"
#include "PicoVetoSteppingAction.hh"
// #include "G4VSensitiveDetector.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " PicoVeto [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "
           << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 9 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif

  G4long myseed = 345354;
  for ( G4int i=1; i<argc; i=i+2 ) {
     if      ( G4String(argv[i]) == "-m" ) macro   = argv[i+1];
     else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
     else if ( G4String(argv[i]) == "-r" ) myseed  = atoi(argv[i+1]);
#ifdef G4MULTITHREADED
     else if ( G4String(argv[i]) == "-t" ) {
                    nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif


// Activate command-based scorer
  // G4ScoringManager::GetScoringManager();
 /////// G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager(); 
 /////// scoringManager->SetVerboseLevel(1);

  
  // Seed the random number generator manually
  G4Random::setTheSeed(myseed);

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager-> SetUserInitialization(new PicoVetoDetectorConstruction());
  // Physics list
  runManager-> SetUserInitialization(new PicoVetoPhysicsList());
  // User action initialization
  runManager->SetUserInitialization(new PicoVetoActionInitialization());

  
 G4UserRunAction* run_action = new PicoVetoRunAction;
  runManager->SetUserAction(run_action);


  
  /////// comment this section cuz :
  // -------- EEEE ------- G4Exception-START -------- EEEE -------
 // *** G4Exception : Run0125
 //       issued by : G4MTRunManager::SetUserAction()
 // For multi-threaded version, define G4UserStackingAction in G4VUserActionInitialization.
 // *** Fatal Exception *** core dump ***
 // -------- EEEE -------- G4Exception-END --------- EEEE -------
 /////// in LINUX but not in wondows ??????? 
  
  
      
  // UserAction classes
  //
  //
  ///  G4UserEventAction* event_action = new PicoVetoEventAction;
  ///  runManager->SetUserAction(event_action);
  //
  ///  G4UserSteppingAction* stepping_action = new PicoVetoSteppingAction;
  ///  runManager->SetUserAction(stepping_action);
  //
  ///  G4UserStackingAction* stacking_action = new PicoVetoStackingAction;
  ///  runManager->SetUserAction(stacking_action);
  //

  

  
  
  // Initialize G4 kernel
  //
  runManager->Initialize();

#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
   
  if ( macro.size() ) {
     // Batch mode
     G4String command = "/control/execute ";
     UImanager->ApplyCommand(command+macro);
  }
    else // Define UI session for interactive mode
  {
#ifdef G4UI_USE
     G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);
#ifdef G4VIS_USE
     UImanager->ApplyCommand("/control/execute vis.mac");
#else
     UImanager->ApplyCommand("/control/execute PicoVeto.in");
#endif
     if (ui->IsGUI())
       UImanager->ApplyCommand("/control/execute gui.mac");
     ui->SessionStart();
     delete ui;
#endif
}

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
