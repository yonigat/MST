/*
 * myDetectorSD.cc
 *
 *  Created on: Mar 27, 2017
 *      Author: adamgeva
 */

#include "myDetectorSD.hh"
#include "detectorHit.hh"
#include "G4SystemOfUnits.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4VUserTrackInformation.hh"
#include "B1TrackInformation.hh"


myDetectorSD::myDetectorSD(G4String name)
: G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
	//history collection name
    G4String HCname = "detectorColl";
    collectionName.insert(HCname);
}


myDetectorSD::~myDetectorSD()
{}


void myDetectorSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new detectorHitsCollection(SensitiveDetectorName,collectionName[0]);
    if (fHCID<0)
    {
    	fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    }
    hce->AddHitsCollection(fHCID,fHitsCollection);
}


G4bool myDetectorSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
	// unmark if required certain amount of energy disposition
    //G4double edep = step->GetTotalEnergyDeposit();
    //if (edep==0.) return true;
	// unmark if required certain particle charge
	//G4double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
	//if (charge==0.) return true;

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    //entering the detector
    if (preStepPoint->GetStepStatus()==fGeomBoundary){

    	G4TouchableHistory* touchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
		G4int copyNo = touchable->GetVolume()->GetCopyNo();
		//G4int ReplicaNum0 = touchable->GetReplicaNumber(0);
		//G4int ReplicaNum1 = touchable->GetReplicaNumber(1);
		//G4int depth = touchable->GetHistoryDepth();

		//G4cout <<"Depth: "<<depth << " ReplicaNum0: "<<ReplicaNum0<<" ReplicaNum1: "<<ReplicaNum1<< " copyNo: "<<copyNo<<G4endl;

		G4double hitTime = preStepPoint->GetGlobalTime();
		G4ThreeVector worldPos = preStepPoint->GetPosition();
		G4ThreeVector localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

		//test: getting information from tracks
		//G4Track* track = step->GetTrack();

		//G4VUserTrackInformation* info = track->GetUserInformation();
		//B1TrackInformation* theInfo = (B1TrackInformation*)info;
		//G4cout << "number of compt: "<< theInfo->GetNumberOfCompton() << " number of Rayl: " << theInfo->GetNumberOfRayl() << G4endl;

		//const G4ParticleDefinition* particle = track->GetParticleDefinition();
		//G4String particleName = particle->GetParticleName();
		//G4double particleMass = particle->GetPDGMass();
		//G4double particleCharge = particle->GetPDGCharge();
		//G4cout << "From Particle definition -- > particle name: "<< particleName << ", Mass: "<< particleMass << ", Charge: " << particleCharge << G4endl;

		//const G4DynamicParticle* dynamicParticle = track->GetDynamicParticle();
		//G4double dynamicParticleTotalEnergy = dynamicParticle->GetTotalEnergy();
		//G4double dynamicParticleKineticEnergy = dynamicParticle->GetKineticEnergy();
		//G4cout << "From Particle dynamicParticle -- > dynamicParticleTotalEnergy: "<< dynamicParticleTotalEnergy << ", dynamicParticleKineticEnergy: "<< dynamicParticleKineticEnergy << G4endl;

		//getting the same information from preStepPoint:
		G4double totalEnergy = preStepPoint->GetTotalEnergy();
		//G4double kineticEnergy = preStepPoint->GetKineticEnergy();
		//G4cout << "From Pre Step Point -- > : Total Energy "<< totalEnergy << ", kineticEnergy: "<< kineticEnergy << G4endl;

	// check if this finger already has a hit
	// this is meant to prevent situation of two hits by the same particle or secondary particles - this code needs to be adjusted to one volume
	/*
	G4int ix = -1;
	for (G4int i=0;i<fHitsCollection->entries();i++)
	{
		if ((*fHitsCollection)[i]->GetID()==copyNo)
		{
			ix = i;
			break;
		}
	}

	if (ix>=0)
		// if it has, then take the earlier time
	{
		if ((*fHitsCollection)[ix]->GetTime()>hitTime)fHHC1ID
		{ (*fHitsCollection)[ix]->SetTime(hitTime); }
	}
	else
		// if not, create a new hit and set it to the collection
	{
		detectorHit* hit = new detectorHit(copyNo,hitTime);
		fHitsCollection->insert(hit);
	}
	*/

		detectorHit* hit = new detectorHit(copyNo);
		hit->SetWorldPos(worldPos);
		hit->SetLocalPos(localPos);
		hit->SetTime(hitTime);
		hit->SetTotalEnergy(totalEnergy);

		fHitsCollection->insert(hit);
    }
    return true;
}
