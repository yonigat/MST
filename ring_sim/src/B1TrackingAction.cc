/*
 * B1TrackingAction.cc
 *
 *  Created on: Apr 23, 2017
 *      Author: adamgeva
 */


#include "B1TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "B1TrackInformation.hh"

B1TrackingAction::B1TrackingAction()
{;}

B1TrackingAction::~B1TrackingAction()
{;}

void B1TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //track just been created - split track already have their parent information copied in the operation
	//could be a new track parentID=0 or photon that was created in flourescent
  if(aTrack->GetUserInformation()==0)
  {
    B1TrackInformation* anInfo = new B1TrackInformation(aTrack);
    G4Track* theTrack = (G4Track*)aTrack;
    theTrack->SetUserInformation(anInfo);
  }
}

//void B1TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
//{
//	//why should we copy the information after the track is killed???
//  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
//  if(secondaries)
//  {
//    B1TrackInformation* info = (B1TrackInformation*)(aTrack->GetUserInformation());
//    size_t nSeco = secondaries->size();
//    if(nSeco>0)
//    {
//      for(size_t i=0;i<nSeco;i++)
//      {
//    	  //todo: could write this method more efficiently
//    	if ((*secondaries)[i]->GetUserInformation()==0)
//    	{
//    		B1TrackInformation* infoNew = new B1TrackInformation(info);
//        	(*secondaries)[i]->SetUserInformation(infoNew);
//    	}
//      }
//    }
//  }
//}


