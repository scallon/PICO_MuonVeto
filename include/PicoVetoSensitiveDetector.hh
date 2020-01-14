
#ifndef PicoVetoSensitiveDetector_h
#define PicoVetoSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"



//class PicoVetoDetectorConstruction;

class PicoVetoSensitiveDetector : public G4VSensitiveDetector
{
 
public:PicoVetoSensitiveDetector(G4String);
	virtual ~PicoVetoSensitiveDetector();

	static int numOfHits;
        static G4int g_N_hits; //global # of hits
        static G4double g_E_total; //total energy deposited into the detector
        static G4double g_E_primary;  //muon's energy 
        static G4String g_E_primary_unit; //unit of muon's energy
        static G4int g_PMT_6_hits[7], g_PMT_3_hits[4], g_PMT_hits; //number of hits in each PMT
        static G4double g_PMT_6_E[7], g_PMT_3_E[4], g_PMT_E;    //energy accumulated in each PMT
 
//	static int numOfEvents;

protected:
	G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
};

#endif // PicoVetoSensitiveDetector_h
