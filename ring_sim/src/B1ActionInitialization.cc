#include "B1ActionInitialization.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1EventAction.hh"
#include "B1RunAction.hh"
#include "B1SteppingAction.hh"
#include "B1TrackingAction.hh"
#include "B1StackingAction.hh"

B1ActionInitialization::B1ActionInitialization(B1DetectorConstruction* geometry)
 : G4VUserActionInitialization(),fdet(geometry)
{}


B1ActionInitialization::~B1ActionInitialization()
{}


void B1ActionInitialization::BuildForMaster() const
{
	SetUserAction(new B1RunAction(fdet));
}


void B1ActionInitialization::Build() const
{
  SetUserAction(new B1PrimaryGeneratorAction);
  SetUserAction(new B1RunAction(fdet));
  //SetUserAction(new B1EventAction);
  SetUserAction(new B1SteppingAction);
  SetUserAction(new B1TrackingAction);
  if (pp->int_map["i_KILL_ELECTRONS"] == 1) {
	  SetUserAction(new B1StackingAction);
  }
}  

