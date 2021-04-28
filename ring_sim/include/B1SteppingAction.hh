/*
 * B1SteppingAction.hh
 *
 *  Created on: Apr 23, 2017
 *      Author: adamgeva
 */

#ifndef B1STEPPINGACTION_HH_
#define B1STEPPINGACTION_HH_

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "ParamsSingleton.hh"


class B1SteppingAction : public G4UserSteppingAction
{
  public:
    B1SteppingAction();
   ~B1SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
  private:
   	ParamsSingleton* pp = ParamsSingleton::Instance();

};


#endif /* B1STEPPINGACTION_HH_ */
