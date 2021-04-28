#ifndef B1EventAction_h
#define B1EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"
#include "ParamsSingleton.hh"


// Event action

class B1EventAction : public G4UserEventAction
{
public:
    B1EventAction();
    virtual ~B1EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

private:
    G4int fdetectorHCID;
   	ParamsSingleton* pp = ParamsSingleton::Instance();

};



#endif
