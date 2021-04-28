/*
 * detectorHit.cc
 *
 *  Created on: Mar 27, 2017
 *      Author: adamgeva
 */

#include "detectorHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<detectorHit>* detectorHitAllocator=0;

detectorHit::detectorHit(G4int i)
: G4VHit(), fId(i), fTime(0.),fLocalPos(0), fWorldPos(0), fTotalEnergy(0)
{}

detectorHit::~detectorHit()
{}

void detectorHit::Print()
{
    G4cout << "  Detector[" << fId << "] : time " << fTime/ns
        << " (nsec) --- local (x,y) " << fLocalPos.x()
        << ", " << fLocalPos.y() << "--- Total Energy: " << fTotalEnergy/keV << G4endl;
}

