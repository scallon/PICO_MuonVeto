
#ifndef PicoVetoEventAction_h
#define PicoVetoEventAction_h 1

#include "G4UserEventAction.hh"
//#include "globals.hh"

class PicoVetoEventAction : public G4UserEventAction
{
  public:
    PicoVetoEventAction();
    virtual ~PicoVetoEventAction();
  

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);

  private:
  //FILE *outputFile;
   
};

#endif


