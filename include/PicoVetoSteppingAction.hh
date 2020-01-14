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
// $Id: PicoVetoSteppingAction.hh 69469 2013-05-05 21:42:35Z ihrivnac $
//
/// \file PicoVetoSteppingAction.hh
/// \brief Definition of the PicoVetoSteppingAction class

#ifndef PicoVetoSteppingAction_h
#define PicoVetoSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

/// Stepping action class
/// 

class PicoVetoSteppingAction : public G4UserSteppingAction
{
  public:
    PicoVetoSteppingAction();
    virtual ~PicoVetoSteppingAction();
    static G4int g_track_ID; 
    static G4int g_N_OP;
    static G4int g_bool_Mu_through_PV; 
    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    G4int fScintillationCounter;
    G4int fCerenkovCounter;
    G4int fEventNumber;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
