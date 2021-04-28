


#ifndef B1BOptrFS_hh
#define B1BOptrFS_hh 1

#include "G4VBiasingOperator.hh"
#include "B1BOptnComptSplitting.hh"
#include "B1BOptnRaylSplitting.hh"
class G4BOptnForceFreeFlight;
class G4BOptnForceCommonTruncatedExp;
class G4VProcess;
class G4BiasingProcessInterface;
class G4ParticleDefinition;
class B1BOptrFSTrackData;
#include <vector>
#include <map>
#include "G4ThreeVector.hh"
#include "ParamsSingleton.hh"


class B1BOptrFS : public G4VBiasingOperator {
public:
  //constructors
  B1BOptrFS(G4String particleToForce,G4String name="FixedSplitting");
  ~B1BOptrFS();

private:
  // -- Mandatory from base class :
  virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess) final;
  virtual G4VBiasingOperation*  ProposeOccurenceBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess) final;
  virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess) final;
  // -- optional methods from base class:
public:
  virtual void           Configure() final;
  virtual void  ConfigureForWorker() final;
  virtual void            StartRun() final;
  virtual void       StartTracking( const G4Track* track ) final;
  virtual void         ExitBiasing( const G4Track*, const G4BiasingProcessInterface* ) final {};
  //virtual void         EndTracking() final;

  // -- operation applied:
//  void OperationApplied( const G4BiasingProcessInterface* callingProcess, G4BiasingAppliedCase biasingCase,
//			 G4VBiasingOperation* operationApplied, const G4VParticleChange* particleChangeProduced ) final;
//  void OperationApplied( const G4BiasingProcessInterface* callingProcess, G4BiasingAppliedCase biasingCase,
//  			 G4VBiasingOperation* occurenceOperationApplied, G4double weightForOccurenceInteraction,
//  			 G4VBiasingOperation* finalStateOperationApplied, const G4VParticleChange* particleChangeProduced ) final;


private:
  G4int                                                                 fFSModelID;
  const G4Track*                                                        fCurrentTrack;
  B1BOptrFSTrackData*                                                   fCurrentTrackData;
  std::map< const G4BiasingProcessInterface*, G4BOptnForceFreeFlight* > fFreeFlightOperations;
  //G4BOptnForceCommonTruncatedExp*                                       fSharedForceInteractionOperation;
  B1BOptnComptSplitting*                                                fComptSplittingOperation;
  B1BOptnRaylSplitting*                                                 fRaylSplittingOperation;
  G4double                                                              fInitialTrackWeight;
  G4bool                                                                fSetup;
  const G4ParticleDefinition*                                           fParticleToBias;

  G4int                          								        fSplittingFactorNp;
  G4int                          								        fSplittingFactorNs;

  ParamsSingleton* pp = ParamsSingleton::Instance();

//  G4bool                        									    fBiasPrimaryOnly;
//  G4bool                         									    fBiasOnlyOnce;
//  G4int 																fBiasNTimes;
//  G4int                 									            fNInteractions;


};

#endif
