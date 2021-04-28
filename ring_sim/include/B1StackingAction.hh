/*
 * B1StackingAction.hh
 *
 *  Created on: Sep 27, 2017
 *      Author: adamgeva
 */

#ifndef B1STACKINGACTION_HH_
#define B1STACKINGACTION_HH_


#include "G4UserStackingAction.hh"
#include "globals.hh"

class EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class B1StackingAction : public G4UserStackingAction
{
  public:
    B1StackingAction();
   ~B1StackingAction();

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

  private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* B1STACKINGACTION_HH_ */
