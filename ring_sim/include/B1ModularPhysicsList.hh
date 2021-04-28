/*
 * B1ModularPhysicsList.hh
 *
 *  Created on: Sep 26, 2017
 *      Author: adamgeva
 */

#ifndef B1MODULARPHYSICSLIST_HH_
#define B1MODULARPHYSICSLIST_HH_


#include "G4VModularPhysicsList.hh"

//sensitive detector
class B1ModularPhysicsList : public G4VModularPhysicsList
{
public:
	B1ModularPhysicsList(G4String);
	virtual ~B1ModularPhysicsList();

    void SetCuts();

private:
    G4double fCutForGamma;

};


#endif /* B1MODULARPHYSICSLIST_HH_ */
