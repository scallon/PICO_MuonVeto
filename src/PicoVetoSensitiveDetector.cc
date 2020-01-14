#include "PicoVetoDetectorConstruction.hh"
#include "PicoVetoSensitiveDetector.hh"
#include "G4VHitsCollection.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4StepPoint.hh"
#include "G4Step.hh"

int PicoVetoSensitiveDetector::numOfHits = 0;
G4int PicoVetoSensitiveDetector::g_N_hits = 0; //zero global # of hits
G4int PicoVetoSensitiveDetector::g_PMT_hits=0;
G4double PicoVetoSensitiveDetector::g_PMT_E=0;
G4double PicoVetoSensitiveDetector::g_E_total=0;
G4double PicoVetoSensitiveDetector::g_E_primary=0;
G4String PicoVetoSensitiveDetector::g_E_primary_unit="   ";

//PicoVetoSensitiveDetector::PicoVetoSensitiveDetector(G4String name, PicoVetoDetectorConstruction *det) :
//	G4VSensitiveDetector(name)
PicoVetoSensitiveDetector::PicoVetoSensitiveDetector(G4String name) :
	G4VSensitiveDetector(name)
{
  //gl.TimeStampSensDet = "Executing " __FILE__ " (compiled " __DATE__ " " __TIME__ ")";
 G4cout<<"-->  Sensitive Detector Added!   " << std::endl;  

  
}
PicoVetoSensitiveDetector::~PicoVetoSensitiveDetector()
{
}

G4bool PicoVetoSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{

       numOfHits++; 
		g_N_hits++;

        G4Track* track = aStep->GetTrack(); //this is track of particle
        G4double globalTime = track->GetGlobalTime(); //time of event
        G4double localTime = track->GetLocalTime(); //time of track
        G4StepPoint* postStepPoint = aStep->GetPostStepPoint(); //this is post step
        G4TouchableHistory* theTouchable_post = (G4TouchableHistory*)(postStepPoint->GetTouchable());
        G4VPhysicalVolume* thePhysical_post = theTouchable_post->GetVolume();
	    G4int copyNo_post=thePhysical_post->GetCopyNo();
        G4String phys_detName_post = thePhysical_post->GetName();

        G4String log_detName_post = thePhysical_post->GetLogicalVolume()->GetName(); //this is name of PMT where hit happened
        G4double depositedEnergy = aStep->GetTotalEnergyDeposit();
        G4double stepLength = aStep->GetStepLength(); 
		G4double photonKinEnergy = track->GetKineticEnergy();
		
         G4StepPoint* postPoint = aStep->GetPostStepPoint();
         G4ThreeVector position = postPoint->GetPosition();
		 
		// G4cout <<"local time : " << localTime << G4endl; 
		 //		 G4cout <<"global time : " << globalTime << G4endl; 

		 
        //g_E_total = g_E_total + depositedEnergy; //total deposited energy 

                      //Hits history 
                      std::ofstream ofs_HitsHistory;
                      G4String HitsHistory_file = "hits.txt";
                      ofs_HitsHistory.open(HitsHistory_file, std::ios::app);
                      //ofs_HitsHistory<<g_N_hits<<" "<<log_detName_post<<"  "<<depositedEnergy*1E6<<"  "<<globalTime/nanosecond<<std::endl;
                      ofs_HitsHistory<<g_N_hits << "\t" << photonKinEnergy*1E6 << " \t" << depositedEnergy*1E6 << " \t " << position.x() << "\t" << position.y() << "\t" << position.z() <<std::endl;
                      ofs_HitsHistory.close();
 
        
    /*    for(int i=1; i<=6; i++)
          {
	    if(PicoVetoDetectorConstruction::g_str_PMT_6_names[i]==log_detName_post)  
              {
                 g_PMT_6_hits[i]++; 
                 g_PMT_6_E[i] = g_PMT_6_E[i]+depositedEnergy;
                 g_E_total = g_E_total + depositedEnergy; //total deposited energy 
                      //SAVING OF THE INFORMATION INTO THE PMTname FILE of energies[eV] and time[ns]
                      
                      std::ofstream ofs_INDIVIDUAL_E_6;
                      G4String INDIVIDUAL_E_6_filename = PicoVetoDetectorConstruction::g_str_PMT_6_names[i] + "_eV_ns_out.dat";
                      ofs_INDIVIDUAL_E_6.open(INDIVIDUAL_E_6_filename, std::ios::app);
                      ofs_INDIVIDUAL_E_6<<g_PMT_6_hits[i]<<"  "<<depositedEnergy*1E6<<"  "<<globalTime/nanosecond<<std::endl;
                      ofs_INDIVIDUAL_E_6.close();
   	      }
	  }*/
	  


        //G4cout<<"Hit in "<<log_detName_post<<std::endl;
        //G4cout<<"deposited E, eV = "<<depositedEnergy*1e6<<std::endl;
        //G4cout<<"last step Length = "<<stepLength/m<<std::endl;

		//G4cout << G4endl<< G4endl<< G4endl<< G4endl;
	//G4cout << "### number of hits SD == " << g_N_hits << G4endl<< G4endl;
       
	return true;
}
