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
/// \file PicoVeto/include/PicoVetoDetectorConstruction.hh
/// \brief Definition of the PicoVetoDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PicoVetoDetectorConstruction_h
#define PicoVetoDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LogicalVolume;
class G4VPhysicalVolume;

class PicoVetoDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    PicoVetoDetectorConstruction();
    virtual ~PicoVetoDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();

  private:
    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

    G4double fTank_x;
    G4double fTank_y;
    G4double fTank_z;

    G4double fBubble_x;
    G4double fBubble_y;
    G4double fBubble_z;




  //Dimensions of Shielding Tank
    G4double tank_inner_R;
    G4double tank_outer_R;
    G4double tank_H;
    G4double tank_start_angle;    
    G4double tank_spanning_angle;
    G4double tank_bottom_inner_R;
    G4double tank_bottom_outer_R;
    G4double tank_bottom_H;
    G4double tank_bottom_start_angle;
    G4double tank_bottom_spanning_angle;
    G4double tank_bottom_Z_shift;
    
    //Dimensions of Liner
    G4double liner_inner_R;
    G4double liner_outer_R;
    G4double liner_H;
    G4double liner_start_angle;    
    G4double liner_spanning_angle;
    G4double liner_top_inner_R;
    G4double liner_top_outer_R;
    G4double liner_top_H;
    G4double liner_top_start_angle;
    G4double liner_top_spanning_angle;
    G4double liner_bottom_inner_R;
    G4double liner_bottom_outer_R;
    G4double liner_bottom_H;
    G4double liner_bottom_start_angle;
    G4double liner_bottom_spanning_angle;
    G4double liner_bottom_Z_shift;
 
    //Dimensions of Water volume
    G4double water_inner_R;
    G4double water_outer_R;
    G4double water_H;
    G4double water_start_angle;
    G4double water_spanning_angle;
    G4double water_level_from_top; 

    //PMT dimensions
    G4double PMT_R;         //radius of PMT photocathode 
    G4double PMT_H;         //thickness of PMT's glass
    G4double PMT_R_sphere;  //radius to mare sphere of the PMT

    //Pressure vessel
    G4double pressure_vessel_inner_R;
    G4double pressure_vessel_outer_R;
    G4double pressure_vessel_H;
    G4double pressure_vessel_start_angle;
    G4double pressure_vessel_spanning_angle;
    G4double pressure_vessel_flange_inner_R;
    G4double pressure_vessel_flange_outer_R;
    G4double pressure_vessel_flange_H;
    G4double pressure_vessel_flange_start_angle;
    G4double pressure_vessel_flange_spanning_angle;
    G4double pressure_vessel_lid_inner_R;
    G4double pressure_vessel_lid_outer_R;
    G4double pressure_vessel_lid_H;
    G4double pressure_vessel_lid_start_angle;
    G4double pressure_vessel_lid_spanning_angle;



 // Helper methods
  void DefineMaterials();
  void SetupGeometry();
  void SetupScoring(G4LogicalVolume* scoringVolume);

  // World logical and physical volumes
  G4LogicalVolume*   expHall_log;
  G4VPhysicalVolume* expHall_phys;
  






  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PicoVetoDetectorConstruction_h*/
