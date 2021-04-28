/*
 * B1SteppingAction.cc
 *
 *  Created on: Apr 23, 2017
 *      Author: adamgeva
 */

#include "B1SteppingAction.hh"
#include "G4VUserTrackInformation.hh"
#include "B1TrackInformation.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include "params.hh"
#include "globalFunctions.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

B1SteppingAction::B1SteppingAction()
:G4UserSteppingAction()
{
	if (pp->int_map["i_VERBOSE_G4NAVIGATOR"]==0){
		//suppressing navigator msgs
		G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetPushVerbosity(0);
	}
}



B1SteppingAction::~B1SteppingAction()
{ }


void B1SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	G4StepPoint* endPoint = aStep->GetPostStepPoint();
	G4StepPoint* startPoint = aStep->GetPreStepPoint();

	G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
	G4Track* theTrack = aStep->GetTrack();
	// G4cout << "Step Number = " << theTrack->GetCurrentStepNumber()  << G4endl;
	G4TouchableHistory* startTouchable = (G4TouchableHistory*)(startPoint->GetTouchable());


	G4String startPhysicalName = startTouchable->GetVolume()->GetName();
	//todo: maybe we need to change the depth
	G4Material* voxel_mat = startPoint->GetMaterial();
//	G4Material* voxel_mat = startTouchable->GetVolume()->GetLogicalVolume()->GetMaterial();

	G4VUserTrackInformation* info = theTrack->GetUserInformation();
	B1TrackInformation* theInfo = (B1TrackInformation*)info;
	// G4cout << "STEPPING" <<  " Mat = " << voxel_mat << " proc name = " << procName << G4endl;
	// aStep->SetLastStepFlag();


	segment currSegmant;
	G4int x = G4RunManager::GetRunManager()->GetNonConstCurrentRun()->GetRunID();

	//we are not counting interaction that occur inside the detector
	if((procName == "compt" || procName == "biasWrapper(compt)") && startPhysicalName!="detectorPixelP") {
		//G4cout<<"We have Compton with prePhysical : " << startPhysicalName <<  G4endl;
		theInfo->AddCompton();

	}
	if((procName == "Rayl" || procName == "biasWrapper(Rayl)") && startPhysicalName!="detectorPixelP") {
		//G4cout<<"We have Rayl with prePhysical : " << startPhysicalName << G4endl;
		theInfo->AddRayl();
	}
	currSegmant.Mat = voxel_mat;
	currSegmant.incidentEnergy = startPoint->GetTotalEnergy();
	currSegmant.scatteredEnergy = endPoint->GetTotalEnergy();
	G4double Len = sqrt(pow((endPoint->GetPosition().x() - startPoint->GetPosition().x()),2) +
			   pow((endPoint->GetPosition().y() - startPoint->GetPosition().y()),2) +
			   pow((endPoint->GetPosition().z() - startPoint->GetPosition().z()),2));
	currSegmant.pathLen = Len/cm; //from mm to cm
	currSegmant.scatteredAngle = angleBetweenVecs(startPoint->GetMomentumDirection(),endPoint->GetMomentumDirection());
	currSegmant.endingProcess = procName;
	
	//TODO: this is only correct for 1 voxel case - need to generalize!
	if (endPoint->GetStepStatus()!=fWorldBoundary){
			if (startPoint->GetPhysicalVolume()->GetName() == "phantom"){
				currSegmant.voxel = startTouchable->GetReplicaNumber(0) + 1;

			}
			else
			{
				currSegmant.voxel = 0;
			}
		}

	//adding segment to list
	theInfo->AddSegment(currSegmant);
	// G4cout << " INC EN: " << startPoint->GetTotalEnergy()/keV << " SCAT EN: " << endPoint->GetTotalEnergy()/keV << " LEN " << Len/cm << "PROC: " << procName  << " VOXEL: " << currSegmant.voxel << G4endl;
	

//	std::cout << "currSegmant.incidentEnergy" << currSegmant.incidentEnergy << " , currSegmant.scatteredEnergy " <<
//			currSegmant.scatteredEnergy << " , currSegmant.pathLen" << currSegmant.pathLen <<
//			" , currSegmant.thetaScatter" << currSegmant.thetaScatter << " , currSegmant.endingProcess" << currSegmant.endingProcess << std::endl;
//


	//is this necessary?
	//theTrack->SetUserInformation(theInfo);

}





