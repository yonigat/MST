#include "B1BOptnComptSplitting.hh"
#include "B1BOptrFS.hh"
#include <vector>
#include "globalFunctions.hh"

#include "G4BiasingProcessInterface.hh"
#include "G4ParticleChangeForGamma.hh"
#include "B1BOptrFSTrackData.hh"
#include "G4VUserTrackInformation.hh"
#include "B1TrackInformation.hh"
#include "G4TrackStatus.hh"
#include "params.hh"
#include <math.h>
#include "G4SystemOfUnits.hh"

//********************
#include "G4EmCalculator.hh"
#include "G4ParticleTable.hh"
//**********************
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1BOptnComptSplitting::B1BOptnComptSplitting(G4String name)
: G4VBiasingOperation(name),
  fSplittingFactorNp(1),
  fSplittingFactorNs(1),
  fParticleChange()
{
  // -- get ID for FS:
  fFSModelID = G4PhysicsModelCatalog::Register("GenBiasFS");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1BOptnComptSplitting::~B1BOptnComptSplitting()
{
}

//G4Track* B1BOptnComptSplitting::GetSplitTrack()
//{
//	if (!fsplitTracksVector.empty()) {
//		G4Track* track = fsplitTracksVector.back();
//		fsplitTracksVector.pop_back();
//		return track;
//	}
//	return nullptr;
//}

G4VParticleChange*
B1BOptnComptSplitting::
ApplyFinalStateBiasing( const G4BiasingProcessInterface* callingProcess,
                        const G4Track*                            track,
                        const G4Step*                              step,
                        G4bool&                                         )
{


  //fetch track data and split factor
  B1BOptrFSTrackData* AuxTrackData = (B1BOptrFSTrackData*)(track->GetAuxiliaryTrackInformation(fFSModelID));
  G4int splitFactor;
  splitFactor = fSplittingFactorNp;

G4cout << "REACHED" << G4endl;
  // -- Collect compt. process (wrapped process) final state:
  G4VParticleChange* processFinalState =
    callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);

  // -- if no splitting requested, let the compton. process to return directly its
  // -- generated final state:
  if ( splitFactor == 1 ) return processFinalState;

  // -- a special case here: the brem. process corrects for cross-section change
  // -- over the step due to energy loss by sometimes "abandoning" the interaction,
  // -- returning an unchanged incoming electron/positron.
  // -- We respect this correction, and if no secondary is produced, its means this
  // -- case is happening:
  //TODO: I think this is relevent for compton process as well. In case the interaction has been cancelled there will not be an electorn as a secondary particle.
  if ( processFinalState->GetNumberOfSecondaries() == 0 )  return processFinalState;

  // -- Now start the biasing:
  // -- We called the compton. process above. Its concrete particle change is indeed
  // -- a "G4ParticleChangeForGamma" object. We cast this particle change to access
  // -- methods of the concrete G4ParticleChangeForGamma type:

  G4ParticleChangeForGamma* actualParticleChange =
    ( G4ParticleChangeForGamma* ) processFinalState ;

  fParticleChange.Initialize(*track);

  //*** -- Store first gamma final state:

  //main photon maintains its weight
  G4double gammaWeight = actualParticleChange->GetWeight();


  fParticleChange.ProposeWeight(gammaWeight);
  //G4cout << "gammaWeight " << gammaWeight << G4endl;
  fParticleChange.ProposeTrackStatus( actualParticleChange->GetTrackStatus() );
  fParticleChange.SetProposedKineticEnergy( actualParticleChange->GetProposedKineticEnergy() );
  fParticleChange.ProposeMomentumDirection( actualParticleChange->GetProposedMomentumDirection() );
  fParticleChange.ProposePolarization( actualParticleChange->GetProposedPolarization() );
  AuxTrackData->fFSState = FSState::start;
  track->SetAuxiliaryTrackInformation(fFSModelID, AuxTrackData);



  // ***-- inform we will have (fSplittingFactor) gamma's + 1 electron:
  fParticleChange.SetNumberOfSecondaries( splitFactor + 1 );

  // -- inform we take care of secondaries weight (otherwise these
  // -- secondaries are by default given the primary weight).
  fParticleChange.SetSecondaryWeightByProcess(true);


  //store electron final state:
   G4Track* electronTrack = actualParticleChange->GetSecondary(0);
   //TODO: is this necessary?
   electronTrack->SetWeight( electronTrack->GetWeight() );
   //G4cout << "electron Weight " << electronTrack->GetWeight() << G4endl;

   //TODO:check that is correct call to mother method?
   fParticleChange.G4VParticleChange::AddSecondary( electronTrack );
   // -- and clean-up the compt. process particle change:
   actualParticleChange->Clear();


  // ***-- now start the fSplittingFactor calls to the compt. process to store each
  // -- related gamma:
  G4int nCalls = 0;
  G4Track* gammaTrack;
  while ( nCalls < splitFactor )
  {
	  //sample again
	  G4VProcess* cal = callingProcess->GetWrappedProcess();
	  G4String nameP = cal->GetProcessName();

	  processFinalState = callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);
	  actualParticleChange = ( G4ParticleChangeForGamma* ) processFinalState ;

      if ( processFinalState->GetNumberOfSecondaries() == 1 )
        {
    	  //TODO: is this OK??? Maybe i need to change also the track information??
    	  gammaTrack = new G4Track( *track );

    	  //TODO: note that the trackInformation did not change! all "original" information is kept, will this affect my scoring system???
    	  //tracking information copy is dealt with in post tracking
    	  //TODO: is the track id being set automatically as a child?
    	  //maybe other values of the tracking information needs to be changed!
//    	  G4VUserTrackInformation* info = track->GetUserInformation();


          gammaTrack->SetWeight( gammaWeight/ splitFactor);
          gammaTrack->SetKineticEnergy(actualParticleChange->GetProposedKineticEnergy());
		  gammaTrack->SetMomentumDirection(actualParticleChange->GetProposedMomentumDirection());
		  gammaTrack->SetTrackStatus(actualParticleChange->GetTrackStatus());
		  gammaTrack->SetPolarization(actualParticleChange->GetProposedPolarization());

		  //check if directed towards to detector
		  G4bool out = outOfRing (gammaTrack->GetPosition(), actualParticleChange->GetProposedMomentumDirection(), pp->double_map["f_DETECTOR_Y"]*mm ,
				  pp->double_map["f_DETECTOR_Y"]*mm*(-1), pp->double_map["f_RADIUS"]*cm);
		  if (!out) //directed at the detector
		  {
			  B1BOptrFSTrackData* SecondaryAuxTrackData = new B1BOptrFSTrackData(AuxTrackData->GetOptr());
			  SecondaryAuxTrackData->fFSState = FSState::toBeFreeFlight;
			  SecondaryAuxTrackData->SetSecondary();
			  gammaTrack->SetAuxiliaryTrackInformation(fFSModelID, SecondaryAuxTrackData);
			  //take regular track information
			  B1TrackInformation* info = (B1TrackInformation*)(track->GetUserInformation());
			  B1TrackInformation* infoNew = new B1TrackInformation(info);
			  infoNew->AddCompton();
			  gammaTrack->SetUserInformation(infoNew);
//			  fsplitTracksVector.push_back(gammaTrack);
			  fParticleChange.G4VParticleChange::AddSecondary( gammaTrack );
          }
		  else//not direct to the detector
	      {
		   //just don't add to secondaries list
			  delete gammaTrack;
	      }
	  }
      // -- very rare special case: we ignore for now.
      //TODO:what is this???
      else if (  processFinalState->GetNumberOfSecondaries() > 1 )
      {
          for ( G4int i = 0 ; i < processFinalState->GetNumberOfSecondaries() ; i++)
            delete processFinalState->GetSecondary(i);
      }

      nCalls++;
      actualParticleChange->Clear();
  }

  // -- we are done:
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
