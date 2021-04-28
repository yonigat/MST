/*
 * B1UserTrackingAction.hh
 *
 *  Created on: Apr 23, 2017
 *      Author: adamgeva
 */

#include "G4UserTrackingAction.hh"

#ifndef B1USERTRACKINGACTION_HH_
#define B1USERTRACKINGACTION_HH_

class B1TrackingAction : public G4UserTrackingAction

{
public:
    B1TrackingAction();
    virtual ~B1TrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
//    virtual void PostUserTrackingAction(const G4Track*);


};





#endif /* B1USERTRACKINGACTION_HH_ */
