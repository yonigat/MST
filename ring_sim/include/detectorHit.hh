/*
 * detectorHit.hh
 *
 *  Created on: Mar 27, 2017
 *      Author: adamgeva
 */

#ifndef DETECTORHIT_HH_
#define DETECTORHIT_HH_

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// detector hit
/// It records:
/// - the strip ID
/// - the particle time
/// - the particle location and energy


class detectorHit : public G4VHit
{
public:
    detectorHit(G4int i);
    virtual ~detectorHit();

    inline void *operator new(size_t);
    inline void operator delete(void*aHit);

    void Print();

    G4int GetID() const { return fId; }

	void SetTime(G4double t) { fTime = t; }
	G4double GetTime() const { return fTime; }

	void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
	G4ThreeVector GetLocalPos() const { return fLocalPos; }

	void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
	G4ThreeVector GetWorldPos() const { return fWorldPos; }

	void SetTotalEnergy(G4double totalEnergy) { fTotalEnergy = totalEnergy; }
	G4double GetTotalEnergy() const { return fTotalEnergy; }

private:
    G4int fId; // this is the copy number
    G4double fTime;
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
    G4double fTotalEnergy;

};

typedef G4THitsCollection<detectorHit> detectorHitsCollection;

extern G4ThreadLocal G4Allocator<detectorHit>* detectorHitAllocator;

inline void* detectorHit::operator new(size_t)
{
    if (!detectorHitAllocator)
        detectorHitAllocator = new G4Allocator<detectorHit>;
    return (void*)detectorHitAllocator->MallocSingle();
}

inline void detectorHit::operator delete(void*aHit)
{
    detectorHitAllocator->FreeSingle((detectorHit*) aHit);
}



#endif /* DETECTORHIT_HH_ */
