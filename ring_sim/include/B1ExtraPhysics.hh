/*
 * B1ExtraPhysics.hh
 *
 *  Created on: Sep 14, 2017
 *      Author: adamgeva
 */

#ifndef B1EXTRAPHYSICS_HH_
#define B1EXTRAPHYSICS_HH_

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

class B1ExtraPhysics : public G4VPhysicsConstructor
{
  public:

    B1ExtraPhysics();
    virtual ~B1ExtraPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

};

#endif


