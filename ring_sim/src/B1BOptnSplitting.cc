#include "B1BOptnSplitting.hh"
#include "B1BOptrFD.hh"
#include <vector>
#include "globalFunctions.hh"

#include "G4BiasingProcessInterface.hh"
#include "G4ParticleChangeForGamma.hh"
#include "B1BOptrFDTrackData.hh"
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

B1BOptnSplitting::B1BOptnSplitting(G4String name)
: G4VBiasingOperation(name),

  fParticleChange()
{
  // -- get ID for FD:
  fFDModelID = G4PhysicsModelCatalog::Register("GenBiasFD");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1BOptnSplitting::~B1BOptnSplitting()
{
}

//G4Track* B1BOptnSplitting::GetSplitTrack()
//{
//	if (!fsplitTracksVector.empty()) {
//		G4Track* track = fsplitTracksVector.back();
//		fsplitTracksVector.pop_back();
//		return track;
//	}
//	return nullptr;
//}

G4VParticleChange*
B1BOptnSplitting::
ApplyFinalStateBiasing( const G4BiasingProcessInterface* callingProcess,
                        const G4Track*                            track,
                        const G4Step*                              step,
                        G4bool&                                         )
{


	//fetch track data and split factor
	B1BOptrFDTrackData* AuxTrackData = (B1BOptrFDTrackData*)(track->GetAuxiliaryTrackInformation(fFDModelID));

	// -- Collect compt. process (wrapped process) final state:
	G4VParticleChange* processFinalState =
	callingProcess->GetWrappedProcess()->PostStepDoIt(*track, *step);

	// -- Now start the biasing:
	// -- We called the compton. process above. Its concrete particle change is indeed
	// -- a "G4ParticleChangeForGamma" object. We cast this particle change to access
	// -- methods of the concrete G4ParticleChangeForGamma type:

	G4ParticleChangeForGamma* actualParticleChange =
	( G4ParticleChangeForGamma* ) processFinalState ;
	G4ThreeVector pos =  track->GetPosition();
	//check if directed towards to detector
	G4bool out = outOfRing (pos, actualParticleChange->GetProposedMomentumDirection(), pp->double_map["f_DETECTOR_Y"]*mm ,
			pp->double_map["f_DETECTOR_Y"]*mm*(-1), pp->double_map["f_RADIUS"]*cm);
	//todo: fixs
	out = false;
	if (out) return processFinalState; //not directed t the detector so return with no splitting

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
	AuxTrackData->fFDState = FDState::normal;
	track->SetAuxiliaryTrackInformation(fFDModelID, AuxTrackData);


	//store secondaries final state:
	G4int number_of_secondaries = actualParticleChange->GetNumberOfSecondaries();

	// -- inform we take care of secondaries weight (otherwise these
	// -- secondaries are by default given the primary weight).
	//fParticleChange.SetSecondaryWeightByProcess(true);

	// ***-- inform we will have one extra secondary
	fParticleChange.SetNumberOfSecondaries( number_of_secondaries + 1 );

	// add other secondaries to particle change
	for (G4int secondary_num = 0; secondary_num < number_of_secondaries; secondary_num++){
	  G4Track* secondary_track = actualParticleChange->GetSecondary(secondary_num);
	  //TODO:check that is correct call to mother method?
	  fParticleChange.G4VParticleChange::AddSecondary( secondary_track );
	}

	//add the extra gamma:
	G4Track* new_gamma_track = new G4Track( *track );

	//TODO: note that the trackInformation did not change! all "original" information is kept, will this affect my scoring system???
	//tracking information copy is dealt with in post tracking
	//TODO: is the track id being set automatically as a child?
	//maybe other values of the tracking information needs to be changed!
	//    	  G4VUserTrackInformation* info = track->GetUserInformation();


	new_gamma_track->SetWeight( gammaWeight );
	new_gamma_track->SetKineticEnergy(actualParticleChange->GetProposedKineticEnergy());
	new_gamma_track->SetMomentumDirection(actualParticleChange->GetProposedMomentumDirection());
	new_gamma_track->SetTrackStatus(actualParticleChange->GetTrackStatus());
	new_gamma_track->SetPolarization(actualParticleChange->GetProposedPolarization());

	B1BOptrFDTrackData* SecondaryAuxTrackData = new B1BOptrFDTrackData(AuxTrackData->GetOptr());
	SecondaryAuxTrackData->fFDState = FDState::toBeFreeFlight;
	SecondaryAuxTrackData->SetSecondary();
	new_gamma_track->SetAuxiliaryTrackInformation(fFDModelID, SecondaryAuxTrackData);

	//take regular track information and update with current interaction
	B1TrackInformation* info = (B1TrackInformation*)(track->GetUserInformation());
	B1TrackInformation* infoNew = new B1TrackInformation(info);

	G4String curr_interaction = callingProcess->GetWrappedProcess()->GetProcessName();
	if (curr_interaction == "compt") infoNew->AddCompton();
	if (curr_interaction == "rayl") infoNew->AddRayl();

	new_gamma_track->SetUserInformation(infoNew);
	fParticleChange.G4VParticleChange::AddSecondary( new_gamma_track );

	actualParticleChange->Clear();
	// -- we are done:
	return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
