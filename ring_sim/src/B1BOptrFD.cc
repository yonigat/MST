
//
#include "B1BOptrFD.hh"
#include "B1BOptrFDTrackData.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4PhysicsModelCatalog.hh"

#include "G4BOptnForceCommonTruncatedExp.hh"
#include "G4ILawCommonTruncatedExp.hh"
#include "G4BOptnForceFreeFlight.hh"
#include "B1BOptnSplitting.hh"


#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4SystemOfUnits.hh"
#include "params.hh"

B1BOptrFD::B1BOptrFD(G4String particleName, G4String name)
  : G4VBiasingOperator(name),
    fFDModelID(-1),
    fCurrentTrack(nullptr),
    fCurrentTrackData(nullptr),
    fInitialTrackWeight(-1.0),
    fSetup(true)
{

  //fSharedForceInteractionOperation = new G4BOptnForceCommonTruncatedExp("SharedForceInteraction");
	fSplittingOperation= new B1BOptnSplitting("Splitting");
  //set splitting factors in these operations
  //todo: is that the best way to do it?

  fParticleToBias = G4ParticleTable::GetParticleTable()->FindParticle(particleName);

  if ( fParticleToBias == 0 )
    {
      G4ExceptionDescription ed;
      ed << " Particle `" << particleName << "' not found !" << G4endl;
      G4Exception(" B1BOptrFD::B1BOptrFD(...)",
		  "BIAS.GEN.07",
		  JustWarning,
		  ed);
    }
}

B1BOptrFD::~B1BOptrFD()
{
  //deleting operations
  for ( std::map< const G4BiasingProcessInterface*, G4BOptnForceFreeFlight* >::iterator it = fFreeFlightOperations.begin() ;
		  it != fFreeFlightOperations.end() ; it++ ) delete (*it).second;
  //delete fSharedForceInteractionOperation;
  delete fSplittingOperation;
}


void B1BOptrFD::Configure()
{
  // -- Create ID for FS:
  fFDModelID = G4PhysicsModelCatalog::Register("GenBiasFD");
  // -- build free flight operations:
  ConfigureForWorker();
}


void B1BOptrFD::ConfigureForWorker()
{
  // -- start by remembering processes under biasing, and create needed biasing operations:
  if ( fSetup )
  {
    // -- get back ID created in master thread in case of MT (or reget just created ID in sequential mode):
    fFDModelID = G4PhysicsModelCatalog::Register("GenBiasFD");

    const G4ProcessManager* processManager = fParticleToBias->GetProcessManager();
    const G4BiasingProcessSharedData* interfaceProcessSharedData = G4BiasingProcessInterface::GetSharedData( processManager );
    if ( interfaceProcessSharedData ) // -- sharedData tested, as is can happen a user attaches an operator
	                                // -- to a volume but without defining BiasingProcessInterface processes.
    {
      for ( size_t i = 0 ; i < (interfaceProcessSharedData->GetPhysicsBiasingProcessInterfaces()).size(); i++ )
	  {
        const G4BiasingProcessInterface* wrapperProcess =
        (interfaceProcessSharedData->GetPhysicsBiasingProcessInterfaces())[i];
        G4String operationName = "FreeFlight-"+wrapperProcess->GetWrappedProcess()->GetProcessName();
	    fFreeFlightOperations[wrapperProcess] = new G4BOptnForceFreeFlight(operationName);
	  }
	}
    fSetup = false;
  }
}


void B1BOptrFD::StartRun()
{
}

G4VBiasingOperation* B1BOptrFD::ProposeNonPhysicsBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{
	G4String process_name = callingProcess->GetProcessName();
	return 0;
}

G4VBiasingOperation* B1BOptrFD::ProposeOccurenceBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{
  // -- does nothing if particle is not of requested type:
  //if ( track->GetDefinition() != fParticleToBias ) return 0;
  //G4int trackID = track->GetTrackID();
  G4String process_name = callingProcess->GetProcessName();

  // -- trying to get auxiliary track data...
//  if ( fCurrentTrackData == nullptr )
//  {
//    // ... and if the track has no aux. track data, it means its biasing is not started yet (note that splitting is to happen first):
//    fCurrentTrackData = (B1BOptrFDTrackData*)(track->GetAuxiliaryTrackInformation(fFSModelID));
//    if ( fCurrentTrackData == nullptr )
//    {
//      return nullptr;
//    }
//  }
  //get aux track data
  fCurrentTrackData = (B1BOptrFDTrackData*)(track->GetAuxiliaryTrackInformation(fFDModelID));
  // -- Send force free flight to the callingProcess:
  // ------------------------------------------------
  // -- The track has been split in the previous step, it has now to be
  // -- forced for a free flight.
  // -- This track will fly with 0.0 weight during its forced flight: ??? maybe not 0??
  // -- it is to forbid the double counting with the force interaction track. no force interaction!
  // -- Its weight is restored at the end of its free flight, this weight
  // -- being its initial weight * the weight for the free flight travel,
  // -- this last one being per process. The initial weight is common, and is
  // -- arbitrary asked to the first operation to take care of it.
  if ( fCurrentTrackData->fFDState == FDState::toBeFreeFlight)
  {
    G4BOptnForceFreeFlight* operation =  fFreeFlightOperations[callingProcess];
    if ( callingProcess->GetWrappedProcess()->GetCurrentInteractionLength() < DBL_MAX/10. )
	{
	  // -- the initial track weight will be restored only by the first DoIt free flight:
	  operation->ResetInitialTrackWeight(track->GetWeight());
      return operation;
    }
  }

  return 0;
}


G4VBiasingOperation* B1BOptrFD::ProposeFinalStateBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{

  //if ( track->GetDefinition() != fParticleToBias ) return 0;
  //fetch track data
  fCurrentTrackData = (B1BOptrFDTrackData*)(track->GetAuxiliaryTrackInformation(fFDModelID));

  //G4int trackID = track->GetTrackID();
  G4String process_name = callingProcess->GetProcessName();

//  if ( fCurrentTrackData != nullptr )
//  {
//     if ( fCurrentTrackData->IsFreeFromBiasing() )
//     {
//  	  // -- takes "ownership" (some track data created before, left free, reuse it). for example if there was Rayl before Compt on the main particle
//	   fCurrentTrackData->fFSOperator = this ;
//     }
//     else
//     {
//	   //in case the particle is free flying
//	   if (fCurrentTrackData->fFSState == FSState::toBeFreeFlight)
//	   {
//		 //!!!!
//         return callingProcess->GetCurrentOccurenceBiasingOperation();
//	   }
//     }
//  }
//  else
//  {
//    fCurrentTrackData = new B1BOptrFDTrackData( this );
//    track->SetAuxiliaryTrackInformation(fFSModelID, fCurrentTrackData);
//  }

  G4VBiasingOperation* operationToReturn = 0;
  //in case the particle is free flying to the detector
  if (fCurrentTrackData->fFDState == FDState::toBeFreeFlight)
  {
	operationToReturn = callingProcess->GetCurrentOccurenceBiasingOperation();
  }
  //in case the particle is free flying and now needs to interact and split
  else if (fCurrentTrackData->fFDState == FDState::normal)
  {
	  //check if its time for compton splitting:
	if (callingProcess->GetProcessName() == "biasWrapper(compt)" || callingProcess->GetProcessName() == "biasWrapper(Rayl)")
	{
	  //send compton splitting operation:
	  // -- Return the compt. splitting operation:
	  fInitialTrackWeight = track->GetWeight();
	  fCurrentTrackData->fFDState = FDState::toBeFreeFlight;
	  operationToReturn = fSplittingOperation;
	}
  }
  //TODO: could the state be any different???

  return operationToReturn;
}


void B1BOptrFD::StartTracking( const G4Track* track )
{
  fCurrentTrack     = track;
//  G4String creating_process = track->GetCreatorModelName();
//  G4String logicalName = track->GetLogicalVolumeAtVertex()->GetName();
//  //G4int ParentID = track->GetParentID();
//  //G4int trackID = track->GetTrackID();
//  G4String volumeName = track->GetVolume()->GetName();
//  G4String nextVolName = track->GetNextVolume()->GetName();
  fCurrentTrackData = (B1BOptrFDTrackData*)track->GetAuxiliaryTrackInformation(fFDModelID);
  if (fCurrentTrackData==nullptr)
  {
	  fCurrentTrackData = new B1BOptrFDTrackData( this );
	  track->SetAuxiliaryTrackInformation(fFDModelID, fCurrentTrackData);
  }
  //for every new track number of compton interaction=0
  //fNInteractions = 0;
}


//void B1BOptrFD::EndTracking()
//{
//}


//void B1BOptrFD::OperationApplied( const G4BiasingProcessInterface*   /*callingProcess*/,
//					      G4BiasingAppliedCase                          BAC,
//					      G4VBiasingOperation*            /* operationApplied*/,
//					      const G4VParticleChange*                          )
//{
//
//  //to be split
//  if  ( fCurrentTrackData->fFSState == FSState::toBeSplitCompt )
//  {
//    fCurrentTrackData->fFSState = FSState::free;
//    G4int numOfCopies = fComptSplittingOperation->GetNumOfTracksCopied();
//    for(G4int i=0;i<numOfCopies;i++)
//    {
//
//      auto cloneData = new B1BOptrFDTrackData( this );
//      cloneData->fFSState = FSState::toBeFreeFlight;
//      //TODO: maybe I should pop the track from the vector to allow another splitting in the future
//   	  fComptSplittingOperation->GetSplitTrack()->SetAuxiliaryTrackInformation(fFSModelID, cloneData);
//    }
//  }
//  else if  ( fCurrentTrackData->fFSState == FSState::toBeSplitRayl )
//    {
//      fCurrentTrackData->fFSState = FSState::free;
//      G4int numOfCopies = fRaylSplittingOperation->GetNumOfTracksCopied();
//      for(G4int i=0;i<numOfCopies;i++)
//      {
//
//        auto cloneData = new B1BOptrFDTrackData( this );
//        cloneData->fFSState = FSState::toBeFreeFlight;
//        //TODO: maybe I should pop the track from the vector to allow another splitting in the future
//     	  fRaylSplittingOperation->GetSplitTrack()->SetAuxiliaryTrackInformation(fFSModelID, cloneData);
//      }
//    }
//  else if ( fCurrentTrackData->fFSState == FSState::toBeFreeFlight )
//  {
//	  //TODO:how can it work with this line?
//      //if ( fFreeFlightOperations[callingProcess]->OperationComplete() ) fCurrentTrackData->Reset(); // -- off biasing for this track
//  }
//  else if ( fCurrentTrackData->fFSState != FSState::start )
//  {
//
//    G4ExceptionDescription ed;
//    ed << " Internal inconsistency : please submit bug report. " << G4endl;
//    G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
//		  "BIAS.GEN.20.4",
//		  JustWarning,
//		  ed);
//  }
//  else if ( fCurrentTrackData->fFSState != FSState::finish )
//  {
//
//      G4ExceptionDescription ed;
//      ed << " Internal inconsistency : please submit bug report. " << G4endl;
//      G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
//  		  "BIAS.GEN.20.4",
//  		  JustWarning,
//  		  ed);
//  }
//  return;
//}
//
//
//void B1BOptrFD::OperationApplied( const G4BiasingProcessInterface* /*callingProcess*/, G4BiasingAppliedCase /*BAC*/,
//			 G4VBiasingOperation* /*occurenceOperationApplied*/, G4double /*weightForOccurenceInteraction*/,
//			 G4VBiasingOperation* /*finalStateOperationApplied*/, const G4VParticleChange* /*particleChangeProduced*/ )
//{
//	  //nothing
//
//}

