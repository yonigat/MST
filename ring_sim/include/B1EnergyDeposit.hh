/*
 * B1EnergyDeposit.hh
 *
 *  Created on: Apr 22, 2017
 *      Author: adamgeva
 */
#include "G4PSEnergyDeposit.hh"
#include "B1TrackInformation.hh"
#include "params.hh"
#include "B1Accumulable.hh"
#include "G4AccumulableManager.hh"
#include "ParamsSingleton.hh"


#ifndef B1ENERGYDEPOSIT_HH_
#define B1ENERGYDEPOSIT_HH_


class B1EnergyDeposit : public G4PSEnergyDeposit
{
  public:
	B1EnergyDeposit(G4String name, G4int type);
    virtual ~B1EnergyDeposit();

    virtual G4int GetIndex(G4Step* step);
	
	G4double weight;
    // grad calculation methods
    G4double getTotalMicXS(G4Element* el, G4double Energy);
    G4double getComptonMicDifferentialXS(G4Element* el, G4double E0 , G4double E1);
    G4double getRaylMicDifferentialXS(G4Element* el, G4double E0, G4double angle);
    G4double getComptonMacDifferentialXS(G4Material* mat, G4double E0 , G4double E1);
    G4double getRaylMacDifferentialXS(G4Material* mat, G4double E0 , G4double angle);
    //update the complete gradient table with current segment contribution
    void updateGradTable(segment seg, G4double final_energy, G4int detIndex);

  protected: // with description
    virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);

  private:
    G4bool recordInteraction (G4Step* aStep,G4TouchableHistory* touchable, G4int totalNumOfInteractions, G4int i, B1TrackInformation* theInfo);
    G4bool recordInteraction_extra (G4Step* aStep,G4TouchableHistory* touchable, G4int totalNumOfInteractions, G4int i, B1TrackInformation* theInfo);
    G4bool recordInteraction_anti_scatter (G4Step* aStep,G4TouchableHistory* touchable, G4double cos_theta);
      //G4int HCID;
      //G4THitsMap<G4double>* EvtMap;
    //scorer type: 0=no_scatter, 1=include_single_scatter, 2=include_multi_scatter, 3=include_single_scatter_compt, 4=include_single_scatter_Rayl
    //5=no_scatter, 6=1_scatter,7=2_scatter,8=3_scatter,9=4_scatter,10=5_scatter,11=6_scatter,12=7_scatter,13=8_scatter,14=9_scatter,15=10_scatter
    //16=with anti-scatter grid!
    G4int fscorerType;
    std::ofstream outputPathsFile;
    // holds the gradient
    B1Accumulable* fGradAccum;
    ParamsSingleton* pp = ParamsSingleton::Instance();

 };



#endif /* B1ENERGYDEPOSIT_HH_ */
