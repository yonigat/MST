/*
 * detectorSD.hh
 *
 *  Created on: Mar 27, 2017
 *      Author: adamgeva
 */

#ifndef DETECTORSD_HH_
#define DETECTORSD_HH_

#include "G4VSensitiveDetector.hh"
#include "detectorHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//sensitive detector
class myDetectorSD : public G4VSensitiveDetector
{
public:
	myDetectorSD(G4String name);
    virtual ~myDetectorSD();
    //HCE = Hits collection of this event
    virtual void Initialize(G4HCofThisEvent* HCE);
    virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist);

private:
    detectorHitsCollection* fHitsCollection;
    //Hits collection ID
    G4int fHCID;
};



#endif /* DETECTORSD_HH_ */
