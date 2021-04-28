

#ifndef B1BOptnRaylSplitting_hh
#define B1BOptnRaylSplitting_hh 1

#include "G4VBiasingOperation.hh"
//#include "G4ParticleChange.hh"
#include "G4ParticleChangeForGamma.hh"
#include <vector>
#include "ParamsSingleton.hh"

class B1BOptnRaylSplitting : public G4VBiasingOperation {
public:
  // -- Constructor :
  B1BOptnRaylSplitting(G4String name);
  // -- destructor:
  virtual ~B1BOptnRaylSplitting();


public:
  // ----------------------------------------------
  // -- Methods from G4VBiasingOperation interface:
  // ----------------------------------------------
  // -- Unused:
  virtual const G4VBiasingInteractionLaw*
  ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*,
                                         G4ForceCondition& )
  { return 0; }

  // --Used:
  virtual G4VParticleChange*   ApplyFinalStateBiasing( const G4BiasingProcessInterface*,
                                                       const G4Track*,
                                                       const G4Step*,
                                                       G4bool&                          );

  // -- Unsued:
  virtual G4double           DistanceToApplyOperation( const G4Track*,
                                                       G4double,
                                                       G4ForceCondition*)
  {return DBL_MAX;}
  virtual G4VParticleChange* GenerateBiasingFinalState( const G4Track*,
                                                        const G4Step*   )
  {return 0;}


public:
  // ----------------------------------------------
  // -- Additional methods, specific to this class:
  // ----------------------------------------------
  // -- Splitting factor:
  void     SetSplittingFactors(G4int splittingFactorNp, G4int splittingFactorNs)
  {
	fSplittingFactorNp = splittingFactorNp;
	fSplittingFactorNs = splittingFactorNs;
  }

  G4int GetSplittingFactorNp() const
  { return fSplittingFactorNp; }
  G4int GetSplittingFactorNs() const
  { return fSplittingFactorNs; }

//  G4Track* GetSplitTrack();
//  G4int GetNumOfTracksCopied()
//  { return fsplitTracksVector.size(); }



private:
  G4int            fSplittingFactorNp;
  G4int            fSplittingFactorNs;
  //G4ParticleChange fParticleChange;
  G4ParticleChangeForGamma fParticleChange;
  //std::vector<G4Track*>  fsplitTracksVector;
  G4int                  fFSModelID;
  ParamsSingleton* pp = ParamsSingleton::Instance();

};

#endif
