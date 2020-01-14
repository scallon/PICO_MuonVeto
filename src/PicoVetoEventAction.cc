//#include <time.h>

//#include <CLHEP/Random/Randomize.h>
//
//#include <stdlib.h>
#define minZ(a,b)  (((a) < (b)) ? (a) : (b)) 

#include "PicoVetoEventAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4VHitsCollection.hh"
#include "G4HCofThisEvent.hh"
#include "G4UserEventAction.hh"

#include "PicoVetoStackingAction.hh"
#include "PicoVetoSteppingAction.hh"
#include "PicoVetoRunAction.hh"
#include "PicoVetoPrimaryGeneratorAction.hh"

#include "PicoVetoSensitiveDetector.hh"

#include "G4Timer.hh"

#include <time.h>     

//#include "LMGlobal.h"
//
//extern LMGlobal gl;

// G4int PicoVetoStackingAction::PMT_detHitsCollectionID = -1; //this global variable declared in StackingAction

PicoVetoEventAction::PicoVetoEventAction()
{//
  //    gl.TimeStampEvAct = "Executing " __FILE__ " (compiled " __DATE__ " " __TIME__ ")";
 
}

PicoVetoEventAction::~PicoVetoEventAction()
{
}

void PicoVetoEventAction::BeginOfEventAction(const G4Event* evt)
{

  G4cout<<std::endl;
  G4cout<<"~~~~~~~~~~~~ Beginning of Event " <<  evt->GetEventID() <<std::endl;

  //G4int event_id = evt->GetEventID();
  //G4cout << "event_id = " << event_id  <<std::endl;
     PicoVetoSensitiveDetector::numOfHits = 0;
     PicoVetoSensitiveDetector::g_N_hits = 0;
     PicoVetoStackingAction::g_N_optical_trajectories = 0;
     //PicoVetoPrimaryGeneratorAction::g_primary_type = "";

     G4SDManager* SDMan = G4SDManager::GetSDMpointer();
     //PicoVetoStackingAction::PMT_detHitsCollectionID = SDMan->GetCollectionID("PMT_detector");
     //G4String sssss = SDMan->GetName(); G4cout<<"----- NNNNNNAMEEEE \n \n \n \n \n \n \n \n \n \n \n \n= "<<sssss<<std::endl;
     //G4cout<<"----------- Hits CollectionID = "<<PicoVetoStackingAction::PMT_detHitsCollectionID<<std::endl;
}

void PicoVetoEventAction::EndOfEventAction(const G4Event* evt)
{

     // timer
//	 G4double timer =  PicoVetoRunAction::myTimer->GetUserElapsed();


  time_t timer;
  G4int seconds;


  time(&timer);  /* get current time; same as: timer = time(NULL)  */
	G4double clock =  PicoVetoRunAction::myTimer; 

	seconds = difftime(timer,clock);
	G4cout << "timer : " <<  (seconds / 60) % 60 << " min " << seconds%60 << " s"<< G4endl; 










 //run #  
       G4int l_run_N = PicoVetoRunAction::g_run_N;
     //  G4cout << "run # = " << l_run_N <<std::endl;

       //event ID = event # in the run
       G4int l_event_ID = evt->GetEventID();
     //  G4cout << "event_ID = " << l_event_ID  <<std::endl;

       //primary particle name !!!!!
       //G4String l_primary_type = "mu-";//PicoVetoPrimaryGeneratorAction::g_primary_type;
       //G4cout << "primary particle : " << l_primary_type << std::endl;

       //# of trajectories
       G4TrajectoryContainer *trajectoryContainer = evt->GetTrajectoryContainer();
       G4int l_N_trajectories = 0;
 //      if (trajectoryContainer) l_N_trajectories = trajectoryContainer->entries();
//       G4cout << "# of trajectories = " << l_N_trajectories<<std::endl; 

       //#ifdef G4VIS_USE
       //for(int i=0; i<=minZ(l_N_trajectories, 100); i++)
       //  {
       //   G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
           //trj->DrawTrajectory();
       // }
       //#endif

       //HITS
       G4HCofThisEvent* PMTHCE = evt->GetHCofThisEvent();
       
       //G4cout<<hc->GetSize()<<"HITS!!!!!!!!!!!!!!!!!!!!"<<std::endl; 

       G4cout<<G4endl; 

       //# of optical photons and trajectories
       G4int l_N_optical_photons = PicoVetoSteppingAction::g_N_OP;  // 
       G4int l_N_optical_trajectories = PicoVetoStackingAction::g_N_optical_trajectories;
 //      G4cout<<"# of optical photons produced until now = "<<l_N_optical_photons<<std::endl;
 //      G4cout<<"# of optical trajectories (optical photons produced in this event) = "<<l_N_optical_trajectories<<std::endl;

       //# of hits
       G4int l_N_hits = PicoVetoSensitiveDetector::g_N_hits;
       G4cout<<"Number of hits on scorer in this event : = "<<l_N_hits<<std::endl;

       //Total energy deposited to the detector, eV
//       G4cout<<"Total enegry deposited to the PMTs = "<<PicoVetoSensitiveDetector::g_E_total*1E6<<" eV"<<std::endl;

       //hits and energies in each PMT
     //  for(int i=1; i<=6; i++)
      //   {
//	   G4cout<<"Number of hits PMT : "<<PicoVetoSensitiveDetector::g_PMT_hits<<std::endl;
//	   G4cout<<"Energy deposited to PMT : "<<PicoVetoSensitiveDetector::g_PMT_E*1E6<<" eV"<<std::endl;
	// }

       G4cout<<G4endl<<G4endl;
/*
       for(int j=1; j<=3; j++)
         {
	   G4cout<<"# of hits PMT_3["<<j<<"] = "<<PicoVetoSensitiveDetector::g_PMT_3_hits[j]<<std::endl;
	   G4cout<<"Energy deposited to PMT_3["<<j<<"] = "<<PicoVetoSensitiveDetector::g_PMT_3_E[j]*1E6<<" eV"<<std::endl;
	 }
*/
/*
       //TOP coordinates of primary particle
       G4double l_rnd_X = PicoVetoPrimaryGeneratorAction::g_rnd_X;
       G4double l_rnd_Y = PicoVetoPrimaryGeneratorAction::g_rnd_Y;
       G4double l_rnd_Z = PicoVetoPrimaryGeneratorAction::g_rnd_Z;
       G4cout<<"Primary TOP coordinates, m:"<<"("<<l_rnd_X<<", "<<l_rnd_Y<<", "<<l_rnd_Z<<")"<<std::endl;

       //Momentum [0, 1]
       G4double l_rnd_X0 = PicoVetoPrimaryGeneratorAction::g_rnd_X0;
       G4double l_rnd_Y0 = PicoVetoPrimaryGeneratorAction::g_rnd_Y0;
       G4double l_rnd_Z0 = PicoVetoPrimaryGeneratorAction::g_rnd_Z0;
       G4cout<<"Primary MOMENTUM: "<<"("<<l_rnd_X0<<", "<<l_rnd_Y0<<", "<<l_rnd_Z0<<")"<<std::endl;
  
       //Angles
       G4double l_theta_deg = PicoVetoPrimaryGeneratorAction::g_theta_deg;
       G4double l_phi_deg   = PicoVetoPrimaryGeneratorAction::g_phi_deg;
       G4cout<<"Primary MOMENTUM angles: "<<"Theta = "<<l_theta_deg<<", "<<"Phi = "<<l_phi_deg<<std::endl;

       //SAVING OF THE INFORMATION INTO THE 'EVENTs' FILE
          std::ofstream ofs_EVENT;
          ofs_EVENT.open("EVENTs_out.dat", std::ios::app);
          ofs_EVENT<<l_run_N<<"   "<<l_event_ID<<"  "<<l_primary_type<<"  "<<l_N_trajectories<<"  "<<l_N_optical_photons<<"  "<<l_N_hits<<std::endl;
          ofs_EVENT.close();

       //SAVING OF THE INFORMATION INTO THE 'COORDINATEs' FILE
          std::ofstream ofs_COORDINATES;
          ofs_COORDINATES.open("COORDINATEs_out.dat", std::ios::app);
          ofs_COORDINATES<<l_run_N<<"  "<<l_event_ID<<"  "<<l_rnd_X<<"  "<<l_rnd_Y<<"  "<<l_rnd_Z<<"  "<<l_rnd_X0<<"  "<<l_rnd_Y0<<"  "<<l_rnd_Z0<<"  "<<l_theta_deg<<"  "<<l_phi_deg<<std::endl;
          ofs_COORDINATES.close();

       //SAVING OF THE INFORMATION INTO THE 'HITs' FILE
          std::ofstream ofs_HIT;
          ofs_HIT.open("HITs_out.dat", std::ios::app);
          ofs_HIT<<l_run_N<<"   "<<l_event_ID<<"  ";
          for(int i=1; i<=6; i++) {ofs_HIT<<PicoVetoSensitiveDetector::g_PMT_6_hits[i]<<"  ";}
          for(int j=1; j<=3; j++) {ofs_HIT<<PicoVetoSensitiveDetector::g_PMT_3_hits[j]<<"  ";}
          ofs_HIT<<l_N_hits<<std::endl;
          ofs_HIT.close();

       //SAVING OF THE INFORMATION INTO THE 'ENERGIEs' FILE [eV]
  
          std::ofstream ofs_ENERGY;
          ofs_ENERGY.open("ENERGIEs[eV]_out.dat", std::ios::app);
          ofs_ENERGY<<l_run_N<<"   "<<l_event_ID<<"  ";
          for(int i=1; i<=6; i++) {ofs_ENERGY<<PicoVetoSensitiveDetector::g_PMT_6_E[i]*1E6<<"  ";}
          for(int j=1; j<=3; j++) {ofs_ENERGY<<PicoVetoSensitiveDetector::g_PMT_3_E[j]*1E6<<"  ";}
          ofs_ENERGY<<PicoVetoSensitiveDetector::g_E_total*1E6<<std::endl;
          ofs_ENERGY.close();
	  
       //SAVING OF THE INFORMATION INTO THE 'muTROUGH_PV' FILE
          if(PicoVetoSteppingAction::g_bool_Mu_through_PV == 1)
	    {
               std::ofstream ofs_muTHROUGH_PV;
               ofs_muTHROUGH_PV.open("muTHROUGH_PV_out.dat", std::ios::app);
               ofs_muTHROUGH_PV<<l_run_N<<"   "<<l_event_ID<<"  ";
               for(int i=1; i<=6; i++) {ofs_muTHROUGH_PV<<PicoVetoSensitiveDetector::g_PMT_6_hits[i]<<"  ";}
               for(int j=1; j<=3; j++) {ofs_muTHROUGH_PV<<PicoVetoSensitiveDetector::g_PMT_3_hits[j]<<"  ";}
               ofs_muTHROUGH_PV<<l_N_hits<<"  "<<l_N_optical_photons<<std::endl;
               ofs_muTHROUGH_PV.close();
	    }

       //SAVING OF THE INFORMATION INTO THE 'PRIMARY particle' FILE
          std::ofstream ofs_PRIMARY;
          ofs_PRIMARY.open("PRIMARY_out.dat", std::ios::app);
          ofs_PRIMARY<<l_run_N<<"  "<<l_event_ID<<"  "<<PicoVetoSensitiveDetector::g_E_primary<<"  "<<PicoVetoSensitiveDetector::g_E_primary_unit<<std::endl;
          ofs_PRIMARY.close();

	  //ZERO GLOBALS
          PicoVetoSensitiveDetector::numOfHits = 0;
          PicoVetoSensitiveDetector::g_N_hits = 0; //zero global # of hits
          PicoVetoSensitiveDetector::g_E_total=0;
          PicoVetoSensitiveDetector::g_E_primary=0;
          PicoVetoSteppingAction::g_N_OP=0; //zero global number of optical photons
          PicoVetoStackingAction::g_N_optical_trajectories = 0; //zero number of optical trajectories
	  PicoVetoSteppingAction::g_bool_Mu_through_PV = 0; //zero flag
           

	  for(int k=0; k<=6; k++) 
            {
               PicoVetoSensitiveDetector::g_PMT_6_hits[k]=0; 
               PicoVetoSensitiveDetector::g_PMT_6_E[k]=0;
            }
          for(int l=0; l<=3; l++) 
            {
               PicoVetoSensitiveDetector::g_PMT_3_hits[l]=0; 
               PicoVetoSensitiveDetector::g_PMT_3_E[l]=0;
            }
*/
 // G4cout<<std::endl;
 // G4cout<<"--------- End of Event ------------"<<std::endl;
 // G4cout<<std::endl;



}
