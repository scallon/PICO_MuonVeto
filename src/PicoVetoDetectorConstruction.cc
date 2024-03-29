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
/// \file PicoVeto/src/PicoVetoDetectorConstruction.cc
/// \brief Implementation of the PicoVetoDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PicoVetoDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4Sphere.hh"
#include "G4VisAttributes.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "PicoVetoSensitiveDetector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoDetectorConstruction::PicoVetoDetectorConstruction()  
 : G4VUserDetectorConstruction() ,expHall_log(0)  ,expHall_phys(0)
{
  //VARIABLES
  fExpHall_x = fExpHall_y = fExpHall_z = 20.0*m;

 
  tank_bottom_H = 0.01*m;                  
  tank_H = (7.9/2)*m; tank_inner_R = 2.8*m; tank_outer_R = 2.85*m;  //Tank cylinder: H=7.9m R=2.8m Thickness=5cm
  tank_start_angle = 0.*deg; tank_spanning_angle = 360.*deg;        //Tank cylinder circle
  tank_bottom_inner_R = 0.0*m; tank_bottom_outer_R = 2.85*m;        //Tank bottom plate: R=2.8m Thickness=5cm  
  tank_bottom_start_angle = 0.*deg; tank_bottom_spanning_angle = 360.*deg;   //Tank bottom circle
  tank_bottom_Z_shift = tank_H+tank_bottom_H;                       //Z position of bottom part of the tank

  water_level_from_top = 30.0*cm;                                //distance from the top of the water tank to the water level
  water_H = tank_H-water_level_from_top/2-(100E-6)*m;               //Water volume
  water_inner_R = 0.0*m;
  water_outer_R = tank_inner_R-(100E-6)*m;  // -0.1 mm to avoid overlapps ??? 
  water_start_angle = 0.*deg;
  water_spanning_angle = 360.*deg; 

  liner_top_inner_R = 0.0*m; 
  liner_top_outer_R = water_outer_R-(100E-6)*m; // as in liner_bottom_outer_R  // (265/2)*cm; // ??
  liner_top_H = 2*cm;
  liner_top_start_angle = 0.*deg;
  liner_top_spanning_angle = 360.*deg;
  //Liner cylinder: H=3.6m-1mm  D=2.8m-1mm R=1.4m-0.5mm 
  liner_bottom_H = (100E-6)*m;
  liner_H = water_H-liner_bottom_H;
  liner_outer_R = water_outer_R-(100E-6)*m;
  liner_inner_R = liner_outer_R-(100E-6)*m;  
  liner_start_angle = 0.*deg;
  liner_spanning_angle = 360.*deg;           //Liner cylinder circle 
  //liner_bottom_H = 0.001*m;
  liner_bottom_inner_R = 0.0*m;
  liner_bottom_outer_R = liner_outer_R;    //Liner botton part: R=1.4m-0.5mm Thickness=1mm
  liner_bottom_start_angle = 0.*deg;
  liner_bottom_spanning_angle = 360.*deg; //Liner bottom part circle
  liner_bottom_Z_shift = liner_H+liner_bottom_H;           //Z position of Liner bottom

  PMT_R = (190/2)*mm;               //PMT radius
  PMT_R_sphere = 131*mm;            //radius of the PMT sphere
  PMT_H = 0.1*mm;                   //PMT thickness 

  
  pressure_vessel_inner_R = (0.77)*m;                      //pressure vessel inner R
  pressure_vessel_outer_R = (0.87)*m;                      //pressure vessel external R
  pressure_vessel_H = (2.4/2)*m;                            //pressure vessel H
  pressure_vessel_start_angle = 0*deg;
  pressure_vessel_spanning_angle = 360*deg;
  
  pressure_vessel_flange_inner_R = pressure_vessel_inner_R;  //pressure vessel flange inner R  !!! at bottom! 
  pressure_vessel_flange_outer_R = (0.90)*m;               //pressure vessel flange external R             
  pressure_vessel_flange_H = (0.20/2)*m;                     //pressure vessel high flange H
  pressure_vessel_flange_start_angle = 0*deg;
  pressure_vessel_flange_spanning_angle = 360*deg;

  pressure_vessel_lid_inner_R = 0*m;                         //pressure vessel lid inner R   !!! at top! 
  pressure_vessel_lid_outer_R = pressure_vessel_outer_R;       // (0.92/2)*m;                  //pressure vessel lid external R             
  pressure_vessel_lid_H = (0.12/2)*m;                        //pressure vessel lid H
  pressure_vessel_lid_start_angle = 0*deg;
  pressure_vessel_lid_spanning_angle = 360*deg;


  

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoDetectorConstruction::~PicoVetoDetectorConstruction(){;}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* PicoVetoDetectorConstruction::Construct()
{
  // Material Definition
  DefineMaterials();  

  // Geometry Definition
  SetupGeometry();   

  // Return world volume
  return expHall_phys;  
}



// ------------------------------------------------------------------------------
// --------------------------------- Materials ----------------------------------
// ------------------------------------------------------------------------------

void PicoVetoDetectorConstruction::DefineMaterials()
{
	
  G4double a, z, density;
  G4int nelements;

//
// Air
//
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  air->AddElement(N, 70.*perCent);
  air->AddElement(O, 30.*perCent);

//
// Water
//
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);


//
// Stainless steel
//
  G4Element* Cr = new G4Element("Chromium"   , "Cr", z=24 , a=52.00*g/mole); 
  G4Element* Mn = new G4Element("Maganese"   , "Mn", z=25 , a=54.94*g/mole); 
  G4Element* Ni = new G4Element("Nickel"     , "Ni", z=28 , a=58.63*g/mole); 
  G4Element* Fe = new G4Element("Iron"       , "Fe", z=26 , a=55.85*g/mole); 
  G4Element* S  = new G4Element("Sulpur"     , "S",  z=16 , a=32.07*g/mole); 
  G4Element* P  = new G4Element("Phosphorus" , "P",  z=15 , a=30.97*g/mole); 
  G4Element* Si = new G4Element("Silicon"    , "Si", z=14 , a=28.08*g/mole); 
  G4Element* C  = new G4Element("Carbon"     , "C",  z=6  , a=12.01*g/mole); 

  G4Material* StainlessSteel302 = new G4Material("StainlessSteel", density=8.03*g/cm3, nelements=8);
  StainlessSteel302->AddElement(Cr,   19.*perCent);
  StainlessSteel302->AddElement(Ni,   10.*perCent);
  StainlessSteel302->AddElement(Mn,    2.*perCent);
  StainlessSteel302->AddElement(C,   0.15*perCent);
  StainlessSteel302->AddElement(S,   0.03*perCent);
  StainlessSteel302->AddElement(P,   0.05*perCent);
  StainlessSteel302->AddElement(Si,    1.*perCent);
  StainlessSteel302->AddElement(Fe, 67.77*perCent);


//
// Polyvinyl chloride PVC flexible
//
  G4Element* Cl = new G4Element("Chlorine"   , "Cr", z=17 , a=35.45*g/mole); 

  G4Material* PVC = new G4Material("PVC", density=1.35*g/cm3, nelements=3);
  PVC->AddElement(C,  2);
  PVC->AddElement(H,  3);
  PVC->AddElement(Cl, 1);

//
// TYVEK 
//
  G4Material* TYVEK = new G4Material("TYVEK", density=1.35*g/cm3, nelements=3);
  TYVEK->AddElement(C,  2);
  TYVEK->AddElement(H,  3);
  TYVEK->AddElement(Cl, 1);

//
// Glass for PMTs
//
  G4Material* Glass = new G4Material("Glass", density=2.2*g/cm3, nelements=2);
  Glass->AddElement(Si, 1);
  Glass->AddElement(O,  2);
 

//
// ------------ Generate & Add Material Properties Table ------------
//
  G4double photonEnergy[] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

//
// Water
//
  G4double refractiveIndex1[] =
            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
              1.346,  1.3465, 1.347,  1.3475, 1.348,
              1.3485, 1.3492, 1.35,   1.3505, 1.351,
              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
              1.3545, 1.355,  1.3555, 1.356,  1.3568,
              1.3572, 1.358,  1.3585, 1.359,  1.3595,
              1.36,   1.3608};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));
/*
  G4double absorption[] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };
*/
		   G4double absorption[] =
           {30.448*m,  40.082*m,  60.329*m,  90.174*m, 120.346*m, 130.889*m,
           150.152*m, 170.241*m, 180.868*m, 200.000*m, 260.316*m, 350.714*m,
           450.455*m, 470.619*m, 520.632*m, 520.632*m, 550.556*m, 520.632*m,
           520.632*m, 470.619*m, 450.455*m, 410.667*m, 370.037*m, 330.333*m,
           300.000*m, 280.500*m, 270.000*m, 240.500*m, 220.000*m, 190.500*m,
           170.500*m, 140.500*m };

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
              7.00, 6.00, 5.00, 4.00 };

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));

  
  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
        ->SetSpline(true);

  //myMPT1->AddProperty("RAYLEIGH",photonEnergy, RayleighWater,     nEntries) ->SetSpline(true);

		
  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);
 
  
  G4double energy_water[] = {
     1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
     1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
     1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
     1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
     1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
     2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
     2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
     2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
     2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
     2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
     3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
     3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
     3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
     4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
     5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  };

  const G4int numentries_water = sizeof(energy_water)/sizeof(G4double);

  //assume 100 times larger than the rayleigh scattering for now.
  G4double mie_water[] = {
     167024.4*m, 158726.7*m, 150742  *m,
     143062.5*m, 135680.2*m, 128587.4*m,
     121776.3*m, 115239.5*m, 108969.5*m,
     102958.8*m, 97200.35*m, 91686.86*m,
     86411.33*m, 81366.79*m, 76546.42*m,
     71943.46*m, 67551.29*m, 63363.36*m,
     59373.25*m, 55574.61*m, 51961.24*m,
     48527.00*m, 45265.87*m, 42171.94*m,
     39239.39*m, 36462.50*m, 33835.68*m,
     31353.41*m, 29010.30*m, 26801.03*m,
     24720.42*m, 22763.36*m, 20924.88*m,
     19200.07*m, 17584.16*m, 16072.45*m,
     14660.38*m, 13343.46*m, 12117.33*m,
     10977.70*m, 9920.416*m, 8941.407*m,
     8036.711*m, 7202.470*m, 6434.927*m,
     5730.429*m, 5085.425*m, 4496.467*m,
     3960.210*m, 3473.413*m, 3032.937*m,
     2635.746*m, 2278.907*m, 1959.588*m,
     1675.064*m, 1422.710*m, 1200.004*m,
     1004.528*m, 833.9666*m, 686.1063*m
  };

  assert(sizeof(mie_water) == sizeof(energy_water));

  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3]={0.99,0.99,0.8};

  myMPT1->AddProperty("MIEHG",energy_water,mie_water,numentries_water)
        ->SetSpline(true);
  myMPT1->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);

  G4cout << "Water G4MaterialPropertiesTable" << G4endl;
  ///  myMPT1->DumpTable();

  water->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator

  water->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  
  
//
// Air
//
  G4double refractiveIndex2[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);

  //G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  //myMPT2->DumpTable();

  air->SetMaterialPropertiesTable(myMPT2);

  
  
}



// ------------------------------------------------------------------------------
// --------------------------------- Geometrie ----------------------------------
// ----------------------------- Volumes and Surfaces ---------------------------
// ------------------------------------------------------------------------------

void PicoVetoDetectorConstruction::SetupGeometry()
{
 
   G4Material* air = G4Material::GetMaterial("Air");
   G4Material* StainlessSteel302 = G4Material::GetMaterial("StainlessSteel");
   G4Material* water = G4Material::GetMaterial("Water");
   G4Material* PVC = G4Material::GetMaterial("PVC");
   G4Material* TYVEK = G4Material::GetMaterial("TYVEK");
  
 
 
 
// 
// The experimental Hall Air
//
  
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  // G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,air,"World",0,0,0);
   expHall_log = new G4LogicalVolume(expHall_box,air,"World",0,0,0);

  //G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
  expHall_phys = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

	
	
	
	
	
	
	
	
// 
// The Stainless Steel Water Tank as a CYLINDER (made from Tube and Solid Tube)
//
  
  //turn the water tank to place it vertical
  G4RotationMatrix* xRot = new G4RotationMatrix;
  xRot->rotateX(90*deg); 
  //set the water tank colour
  //G4Colour waterTank_colour(.5,.5,.5); //GREY
  G4Colour waterTank_colour(1.0, 0.0, 1.0); //Magenta
  G4VisAttributes* waterTank_graph = new G4VisAttributes(waterTank_colour);
  waterTank_graph->SetVisibility(true); //colour of Tank is Blue
  //waterTank_graph->SetForceSolid(true); //Visualiszation method of Tank is Solid
  waterTank_graph->SetForceWireframe(true); //Visualiszation method of Tank is Solid
  
  G4Tubs* waterTank_cylinder
    = new G4Tubs("Tank_cylinder", tank_inner_R, tank_outer_R, tank_H, tank_start_angle, tank_spanning_angle);
  G4VSolid* waterTank_bottom 
    = new G4Tubs("Tank_bottom", tank_bottom_inner_R, tank_bottom_outer_R, tank_bottom_H, tank_bottom_start_angle, tank_bottom_spanning_angle);

  G4UnionSolid* waterTank_box
    = new G4UnionSolid("Tank", waterTank_cylinder, waterTank_bottom, 0, G4ThreeVector(0,0,-tank_bottom_Z_shift));
  
  G4LogicalVolume* waterTank_log
    = new G4LogicalVolume(waterTank_box, StainlessSteel302, "Tank",0,0,0);
  waterTank_log->SetVisAttributes(waterTank_graph); 

  G4VPhysicalVolume* waterTank_phys
    = new G4PVPlacement(xRot,G4ThreeVector(0, water_level_from_top/2, 0),waterTank_log,"Tank",
			expHall_log,false,0);

	


//
// The Water volume
//
  
  G4Colour water_colour(0., 0., 1.); //Blue
  G4VisAttributes* water_graph = new G4VisAttributes(water_colour);  
  water_graph->SetVisibility(true); //colour of Water is blue
  //water_graph->SetForceSolid(true); //Visualiszation method of Water is Solid
  water_graph->SetForceWireframe(true); //Visualiszation method of Water is Solid

  G4VSolid* Water_volume 
    = new G4Tubs("Water_volume", water_inner_R, water_outer_R, water_H, water_start_angle, water_spanning_angle);

  G4LogicalVolume* Water_log
    = new G4LogicalVolume(Water_volume, water, "Water",0,0,0);
  Water_log->SetVisAttributes(water_graph); 

   G4VPhysicalVolume* Water_phys
    = new G4PVPlacement(xRot,G4ThreeVector(),Water_log,"Water",
			expHall_log,false,0);




//
// The PVC Liner -------------- WE CAN DECLARE LINER HERE OUT of Water
//
  //set the PVC liner colour
  G4Colour linerPVC_colour(1.,1.,1.); //WHITE
  G4VisAttributes* linerPVC_graph = new G4VisAttributes(linerPVC_colour);  
  //linerPVC_graph->SetVisibility(true); //colour of PVC liner is White
  //linerPVC_graph->SetForceSolid(true); //Visualiszation method of PVC liner is Solid
  linerPVC_graph->SetForceWireframe(true); //Visualiszation method of PVC liner is Solid


  G4Tubs* linerPVC_cylinder
    = new G4Tubs("linerPVC_cylinder", liner_inner_R, liner_outer_R, liner_H, liner_start_angle, liner_spanning_angle);

  G4VSolid* linerPVC_bottom 
    = new G4Tubs("linerPVC_bottom", liner_bottom_inner_R, liner_bottom_outer_R, liner_bottom_H, liner_bottom_start_angle, liner_bottom_spanning_angle);

  
   
  G4UnionSolid* linerPVC_box
    = new G4UnionSolid("Liner", linerPVC_cylinder, linerPVC_bottom, 0, G4ThreeVector(0,0,(-liner_bottom_Z_shift)));

  G4LogicalVolume* linerPVC_log
    = new G4LogicalVolume(linerPVC_box, PVC, "Liner",0,0,0);
  linerPVC_log->SetVisAttributes(linerPVC_graph); 

  G4VPhysicalVolume* linerPVC_phys
       = new G4PVPlacement(0,G4ThreeVector(),linerPVC_log,"Liner",
     			Water_log,false,0);
 











// STEEL Pressure vessel
  // Flange at bottom!  
  // Lid on top! 
  
  //set the pressure vessel colour
  G4Colour pressure_vessel_colour(0.5,0.5,0.5); //GRAY
  G4VisAttributes* pressure_vessel_graph = new G4VisAttributes(pressure_vessel_colour);  
  pressure_vessel_graph->SetVisibility(true); 
  pressure_vessel_graph->SetForceSolid(true); 
  //pressure_vessel_graph->SetForceWireframe(true); 


  // G4Tubs("Tracker", innerRadius, outerRadius, halfLength, startAngle, spanningAngle);
  G4Tubs* pressureVessel_cylinder
    = new G4Tubs("Pressure_vessel_cylinder", pressure_vessel_inner_R, pressure_vessel_outer_R, pressure_vessel_H, pressure_vessel_start_angle, pressure_vessel_spanning_angle);
  G4VSolid* pressureVessel_bottom 
    = new G4Tubs("Pressure_vessel_bottom", 0, pressure_vessel_inner_R, 0.03*m, pressure_vessel_start_angle, pressure_vessel_spanning_angle); 
  G4VSolid* pressureVessel_flange 
    = new G4Tubs("Pressure_vessel_flange", pressure_vessel_flange_inner_R, pressure_vessel_flange_outer_R, pressure_vessel_flange_H, pressure_vessel_flange_start_angle, pressure_vessel_flange_spanning_angle); 
  G4VSolid* pressureVessel_lid 
    //   = new G4Tubs("Pressure_vessel_lid", pressure_vessel_lid_inner_R, pressure_vessel_lid_outer_R, pressure_vessel_lid_H, pressure_vessel_lid_start_angle, pressure_vessel_lid_spanning_angle); 
    //  G4Sphere(const G4String& pName, G4double pRmin, G4double pRmax, G4double pSPhi, G4double pDPhi, G4double pSTheta, G4double pDTheta );
    //  G4Ellipsoid(const G4String& pName, G4double pxSemiAxis, G4double pySemiAxis, G4double pzSemiAxis, G4double pzBottomCut, G4double pzTopCut)
    = new G4Ellipsoid("Pressure_vessel_lid", pressure_vessel_outer_R, pressure_vessel_outer_R, 0.40*m, 0.*m, 0.30*m);

  //  = new G4Sphere("Pressure_vessel_lid", pressure_vessel_inner_R, pressure_vessel_outer_R, 0.*degree, 360.*degree, 0.*degree, 90.*degree);


    
  //  G4UnionSolid (const G4String &pName, G4VSolid *pSolidA, G4VSolid *pSolidB, G4RotationMatrix *rotMatrix, const G4ThreeVector &transVector)
  G4UnionSolid* pressureVessel_cylinder_bottom
    = new G4UnionSolid("pressure_vessel_cylinder_bottom", pressureVessel_cylinder, pressureVessel_bottom, 0, G4ThreeVector(0,0,-(pressure_vessel_H-0.03*m) ));
  G4UnionSolid* pressureVessel_cylinder_bottom_flange
    //    = new G4UnionSolid("pressure_vessel_cylinder_bottom_flange", pressureVessel_cylinder_bottom, pressureVessel_flange, 0, G4ThreeVector(0,0, (pressure_vessel_H+pressure_vessel_flange_H) ));
   = new G4UnionSolid("pressure_vessel_cylinder_bottom_flange", pressureVessel_cylinder_bottom, pressureVessel_flange, 0, G4ThreeVector(0,0, -(pressure_vessel_H+pressure_vessel_flange_H) ));
  G4UnionSolid* pressureVessel_box  
    //    = new G4UnionSolid("pressure_vessel_box", pressureVessel_cylinder_bottom_flange, pressureVessel_lid, 0, G4ThreeVector(0,0, (pressure_vessel_H+2*pressure_vessel_flange_H) ));
   = new G4UnionSolid("pressure_vessel_box", pressureVessel_cylinder_bottom_flange, pressureVessel_lid, 0, G4ThreeVector(0,0, pressure_vessel_H));


  G4LogicalVolume* pressureVessel_log
    = new G4LogicalVolume(pressureVessel_box, StainlessSteel302, "Pressure_vessel",0,0,0);
  pressureVessel_log->SetVisAttributes(pressure_vessel_graph); 


  /*
  // warning: unused variable pressureVessel_phys [-Wunused-variable]
  G4VPhysicalVolume* pressureVessel_phys
    = new G4PVPlacement(0,G4ThreeVector(),pressureVessel_log,"Pressure_vessel",
			Water_log,false,0);
  
  */


  // Legs
  /*
G4CutTubs( const G4String& pName,
                    G4double pRMin,
                    G4double pRMax,
                    G4double pDz,
                    G4double pSPhi,
                    G4double pDPhi,
                    G4ThreeVector pLowNorm,
                    G4ThreeVector pHighNorm )
    

 G4CutTubs* pv_leg_cylinder
  = new G4CutTubs("pv_legs", 0, pv_leg_R, pv_leg_length, pv_leg_start_angle, pv_leg_spanning_angle, G4ThreeVector(),  G4ThreeVector()  );

G4LogicalVolume* pc_leg_log1
    = new G4LogicalVolume(pv_leg_cylinder, StainlessSteel302, "Pressure_vessel_leg1" , 0, 60.96 , 0);
G4LogicalVolume* pc_leg_log2
    = new G4LogicalVolume(pv_leg_cylinder, StainlessSteel302, "Pressure_vessel_leg2" , 0,  , 0);
G4LogicalVolume* pc_leg_log3
    = new G4LogicalVolume(pv_leg_cylinder, StainlessSteel302, "Pressure_vessel_leg3" , 0,  , 0);

  */




  












  // G4UnionSolid (const G4String &pName, G4VSolid *pSolidA, G4VSolid *pSolidB, G4RotationMatrix *rotMatrix, const G4ThreeVector &transVector)

  
// Cylindrical TYVEK around the pressure vessel !!!
  //set the pressure vessel color
  //	  G4Colour TYVEK_colour(1.,1.,0.); //YELLOW
	   G4Colour TYVEK_colour(0.,1.,1.); //CYAN
 G4VisAttributes* TYVEK_graph = new G4VisAttributes(TYVEK_colour);  
  TYVEK_graph->SetVisibility(true); 
  // TYVEK_graph->SetForceSolid(true); //Visualiszation method of TYVEK is Solid
  TYVEK_graph->SetForceWireframe(true); //Visualiszation method of TYVEK is Solid
  
  G4double TYVEK_thickness = 10*mm;
  G4double TYVEK_inner_R = pressure_vessel_flange_outer_R+TYVEK_thickness;
  G4double TYVEK_outer_R = TYVEK_inner_R+TYVEK_thickness;
  G4double TYVEK_H = pressure_vessel_H+pressure_vessel_flange_H+pressure_vessel_lid_H+1.0*mm+30.0*cm; // 30 cm is the pv top thickness 

  G4Tubs* TYVEK_cylinder
    = new G4Tubs("TYVEK_cylinder", TYVEK_inner_R, TYVEK_outer_R, TYVEK_H, liner_start_angle, liner_spanning_angle);
   G4VSolid* TYVEK_bottom 
    = new G4Tubs("TYVEK_bottom", 0, TYVEK_outer_R, TYVEK_thickness, 0*deg, 360*deg);
  G4VSolid* TYVEK_top 
    = new G4Tubs("TYVEK_top", 0, TYVEK_outer_R, TYVEK_thickness, 0*deg, 360*deg);
  
  G4UnionSolid* TYVEK_cylinder_bottom
    = new G4UnionSolid("TYVEK_cylinder_bottom", TYVEK_cylinder, TYVEK_bottom, 0, G4ThreeVector(0,0, -TYVEK_H-TYVEK_thickness));
   G4UnionSolid* TYVEK_box
    = new G4UnionSolid("TYVEK_box", TYVEK_cylinder_bottom, TYVEK_top, 0, G4ThreeVector(0,0, TYVEK_H+TYVEK_thickness));
  
  G4LogicalVolume* TYVEK_log
      = new G4LogicalVolume(TYVEK_box, TYVEK, "TYVEK_LOG",0,0,0);

  TYVEK_log->SetVisAttributes(TYVEK_graph); 

  /*
  //  warning: unused variable TYVEK_phys [-Wunused-variable]
  G4VPhysicalVolume* TYVEK_phys
    = new G4PVPlacement(0,G4ThreeVector(),TYVEK_log,"TYVEK_phys",
			Water_log,false,0);




  
  
          //ATTENTION ! This is the AIR cylinder inside of the TYVEK cylinder to avoid an optical photons inside the cylinder
             
         //  G4Colour TYVEK_302_colour(0.,1.,1.); //CYAN
	  G4Colour TYVEK_302_colour(1.,1.,0.); //YELLOW
          G4VisAttributes* TYVEK_302_graph = new G4VisAttributes(TYVEK_302_colour);  
          TYVEK_302_graph->SetVisibility(true); 
          //TYVEK_302_graph->SetForceSolid(true); //Visualiszation method of TYVEK is Solid
	   TYVEK_302_graph->SetForceWireframe(true); //Visualiszation method of TYVEK is wire

	  
          G4Tubs* TYVEK_302_cylinder  //Steel-302 cylinder to place it in the TYVEK cylinder
                  = new G4Tubs("TYVEK_302_cylinder", 0, TYVEK_inner_R, TYVEK_H-0.01*m, 0*deg, 360*deg);  		       
          G4LogicalVolume* TYVEK_302_log
                  = new G4LogicalVolume(TYVEK_302_cylinder, air, "TYVEK_302",0,0,0);
          TYVEK_302_log->SetVisAttributes(TYVEK_302_graph);
	  
	  // warning: unused variable TYVEK_302_phys [-Wunused-variable]
	  G4VPhysicalVolume* TYVEK_302_phys
                  = new G4PVPlacement(0,G4ThreeVector(),TYVEK_302_log,"TYVEK_302_phys",
                                      Water_log,false,0);  
  	  
  */

  

	  
	  
	  
	  
	  
	  
	
	
//
// Top Box Scorer On Water Tank : 
//
  G4ThreeVector positionTopScorer = G4ThreeVector(0.*cm, 400.*cm, 0.*cm);    
  
  G4Box* solidTopScorer = new G4Box("TopScorerBox", 3.0*m, 0.001*m, 3.0*m); 
 
  G4LogicalVolume* logicTopScorer =  new G4LogicalVolume(solidTopScorer, StainlessSteel302, "TopScorer", 0,0,0);
              
  G4VPhysicalVolume* TopScorer_phys =  new G4PVPlacement(0,positionTopScorer,logicTopScorer,"TopScorer",expHall_log,false,0);  

  //  G4VisAttributes* scoreAttributes =  new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.4));
  G4VisAttributes* scoreAttributes =  new G4VisAttributes(G4Colour(0.0,0.0,0.0,0.)); // black for no show
   scoreAttributes->SetVisibility(true);
   logicTopScorer->SetVisAttributes(scoreAttributes);
  
  
  
  
 //
 // Choose between top scorer or Sides & bottom scorer : 
 //
   // setup scoring Top  
   SetupScoring(logicTopScorer);
   // setup scoring Sides   
  ///  SetupScoring(linerPVC_log);
 
	  
	  
	  
	  
	  


  

// ------------------------------------------------------------------------------
// --------------------------------- Surfaces -----------------------------------
// ------------------------------------------------------------------------------
	  
//
/// Water Props :
//
	  const G4int num = 5; //number for Optical properties
	  G4double ephoton[num] = {2.034*eV, 2.75*eV, 2.5*eV, 3.181*eV, 4.136*eV};

	  //Optical Water Properties
       	  G4double refractiveIndex[num] = {1.3435, 1.3522, 1.35, 1.3555, 1.3608};
	  //G4double refractiveIndex[num] = {1., 1., 1., 1., 1.};
	  G4double specularLobe[num]    = {0.3, 0.3, 0.3, 0.3, 0.3};
	  G4double specularSpike[num]   = {0.2, 0.2, 0.2, 0.2, 0.2};
	  //  G4double specularSpike[num]   = {0.99, 0.99, 0.99, 0.99, 0.99};
	  G4double backScatter[num]     = {0.2, 0.2, 0.2, 0.2, 0.2};
	  
	  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
	  G4MaterialPropertiesTable* myST2 = new G4MaterialPropertiesTable();

	  myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num);

	  myST2->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
	  myST2->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
	  myST2->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
	  myST2->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);


	  
//
// Water Top Surface Border with air
//
  
  G4OpticalSurface* opWaterTopSurface = new G4OpticalSurface("WaterTopSurface");
  opWaterTopSurface->SetType(dielectric_dielectric);
  opWaterTopSurface->SetFinish(polished);
  opWaterTopSurface->SetModel(unified);

  //  new G4LogicalBorderSurface("WaterSurface",waterTank_phys,expHall_phys,opWaterSurface);
new G4LogicalBorderSurface("WaterTopSurface",Water_phys,expHall_phys,opWaterTopSurface);
   G4double ReflectivityAir[num] = {0.01, 0.01, 0.01, 0.01, 0.01}; 
    myST1->AddProperty("REFLECTIVITY", ephoton, ReflectivityAir, num);
    opWaterTopSurface -> SetMaterialPropertiesTable(myST1);





    
//
// Water - PVC (bottom and sides) Border between water and PVC
//
 G4OpticalSurface* OpWaterPVCSurface = new G4OpticalSurface("WaterPVCSurface");
  OpWaterPVCSurface->SetType(dielectric_dielectric);   //two dielectrics water and PVC
  OpWaterPVCSurface->SetFinish(polishedfrontpainted);  //this is type of surface
 //   OpWaterPVCSurface->SetType(dielectric_metal);   // dielectrics water and "metal" PVC
 //  OpWaterPVCSurface->SetFinish(polished);  //this is type of surface
 // OpWaterPVCSurface->SetModel(glisur);                 //model type
  OpWaterPVCSurface->SetModel(unified);                 //model type

  
  //OpWaterPVCSurface -> SetMaterialPropertiesTable(myST2);

  
  //*---------- DEAP reflectivity, efficiency ---------*
  
  //  G4double Reflectivity[num] = {0.81, 0.78, 0.76, 0.2, 0.28}; //Reflectivity depending on wavelength DEAP
  // G4double Efficiency[num]   = {0.9, 0.9, 0.9, 0.9, 0.9};
  
  //---------- Manual reflectivity, efficiency -------
   G4double Reflectivity[num] = {0., 0., 0., 0., 0.}; 
 // G4double Reflectivity[num] = {0.1, 0.1, 0.1, 0.1, 0.1}; 
// G4double Reflectivity[num] = {0.5, 0.5, 0.5, 0.5, 0.5}; 
  // G4double Reflectivity[num] = {0.99, 0.99, 0.99, 0.99, 0.99}; 
  // G4double Reflectivity[num] = {0.999, 0.999, 0.999, 0.999, 0.999}; 

  ///// G4double Reflectivity[num] = {1., 1., 1., 1., 1.}; 

   G4double Efficiency[num]  = {0.9, 0.9, 0.9, 0.9, 0.9};
  //    G4double Efficiency[num]   = {0.01, 0.01, 0.01, 0.01, 0.01};
  //  G4double Efficiency[num]   = {1.0, 1.0, 1.0, 1.0, 1.0};
 // G4double Efficiency[num]   = {0., 0., 0., 0., 0.};


    // G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

   myST1->AddProperty("REFLECTIVITY", ephoton, Reflectivity, num);
    myST1->AddProperty("EFFICIENCY",   ephoton, Efficiency,   num);
  
    OpWaterPVCSurface->SetMaterialPropertiesTable(myST1);

 new G4LogicalBorderSurface("WaterPVCSurface",Water_phys,linerPVC_phys,OpWaterPVCSurface);

 G4cout << "PVC Surface G4MaterialPropertiesTable" << G4endl;
 //// myST1->DumpTable();


 /*


//
// TYVEK - water Border
//
  
  G4OpticalSurface* OpTYVEKSurface = new G4OpticalSurface("TYVEKSurface");
  OpTYVEKSurface->SetType(dielectric_dielectric);   //two dielectrics water and tyvek
  OpTYVEKSurface->SetFinish(polishedfrontpainted);  //this is type of surface
  OpTYVEKSurface->SetModel(unified);                 //model type

  G4LogicalSkinSurface* TYVEKSurface = 
          new G4LogicalSkinSurface("TYVEKSurface", TYVEK_log, OpTYVEKSurface);

  // warning: unused variable opticalTYVEKSurface [-Wunused-variable]
  G4OpticalSurface* opticalTYVEKSurface = dynamic_cast <G4OpticalSurface*>
    (TYVEKSurface->GetSurface(TYVEK_log)->GetSurfaceProperty());

 
  //*------- TYVEK reflectivity, efficiency --------------*
  
  // G4double tyvekReflectivity[num] = {0.9, 0.9, 0.9, 0.885, 0.75}; //Reflectivity depending on wavelength TYVEK
  /////       G4double tyvekReflectivity[num]   = {1.0, 1.0, 1.0, 1.0, 1.0};
   G4double tyvekReflectivity[num]   = {0., 0., 0., 0., 0.};
  G4double tyvekEfficiency[num]   = {0.9, 0.9, 0.9, 0.9, 0.9};
 // G4double tyvekEfficiency[num]   = {0., 0., 0., 0., 0.};
  
 

  G4MaterialPropertiesTable *mySTtyvek = new G4MaterialPropertiesTable();

  mySTtyvek->AddProperty("REFLECTIVITY", ephoton, tyvekReflectivity, num);
  mySTtyvek->AddProperty("EFFICIENCY",   ephoton, tyvekEfficiency,   num);

  OpTYVEKSurface->SetMaterialPropertiesTable(mySTtyvek);

  //  new G4LogicalSkinSurface("TYVEKsurface", TYVEK_log, OpTYVEKSurface);
  // new G4LogicalBorderSurface("TYVEKsurface", TYVEK_phys,Water_phys, OpTYVEKSurface);
new G4LogicalBorderSurface("TYVEKsurface", Water_phys, TYVEK_phys, OpTYVEKSurface);

 */

// Pressure vessel !!! 
//
  // warning: unused variable LinerWater [-Wunused-variable]
  //   G4LogicalBorderSurface* LinerWater = new G4LogicalBorderSurface("LinerWater", linerPVC_phys, Water_phys, opWaterSurface);           //between wall liner, bottom liner and water


  // G4LogicalBorderSurface* TopLinerWater = new G4LogicalBorderSurface("TopLinerWater", linerPVC_top_phys, Water_phys, opWaterSurface); //between top liner and water

  //  G4LogicalBorderSurface::DumpInfo();





//
// Top Scorer
// 

 G4OpticalSurface* OpScorerSurface = new G4OpticalSurface("ScorerSurface");
  OpScorerSurface->SetType(dielectric_metal);   //two dielectrics water and PVC
  OpScorerSurface->SetFinish(polished);  //this is type of surface
  OpScorerSurface->SetModel(unified);                 //model type

  
  G4double ReflectivityScorer[num] = {0., 0., 0., 0., 0.}; 
//  G4double EfficiencyScorer[num]  = {0.025, 0.18, 0.20, 0.25, 0.08}; //this is real PMT efficiencies
   G4double EfficiencyScorer[num]  = {1., 1., 1., 1., 1.};

   G4MaterialPropertiesTable *mySScorer = new G4MaterialPropertiesTable();

   mySScorer->AddProperty("REFLECTIVITY", ephoton, ReflectivityScorer, num);
    mySScorer->AddProperty("EFFICIENCY",   ephoton, EfficiencyScorer,   num);

    OpScorerSurface->SetMaterialPropertiesTable(mySScorer);

 new G4LogicalSkinSurface("ScorerSurface",logicTopScorer,OpScorerSurface);







 

}







// ------------------------------------------------------------------------------
// -------------------------- Sensitive Detector --------------------------------
// ------------------------------------------------------------------------------
 
void PicoVetoDetectorConstruction::SetupScoring(G4LogicalVolume* scoringVolume)
{

  // Get pointer to detector manager
  G4SDManager* manager = G4SDManager::GetSDMpointer();  
 
 
 
  PicoVetoSensitiveDetector* sensitiveDetector = new PicoVetoSensitiveDetector("PMT_detector");

     // Register detector with manager

   manager->AddNewDetector(sensitiveDetector);

   
  // Attach detector to scoring volume
  scoringVolume->SetSensitiveDetector(sensitiveDetector);
 
  
 
 
 
 
 //  G4cout<<"+++++++++++++++++++++++++ SensitiveDetector ----------------------------"<<G4endl<<G4endl<<G4endl<<G4endl<<G4endl;

 
 
 
 
 
 
 
 
 
 
  
}





