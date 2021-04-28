/*
 * B1Run.cc
 *
 *  Created on: Apr 22, 2017
 *      Author: adamgeva
 */


#include "B1Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "globalFunctions.hh"


B1Run::B1Run()
 : G4Run()
{
	fMapSum = new G4THitsMap<G4double> [pp->int_map["i_NUM_OF_SCORERS"]];
	fColIDSum = new G4int [pp->int_map["i_NUM_OF_SCORERS"]];
	G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	G4String fullName;
	for(G4int i=0; i<(pp->int_map["i_NUM_OF_SCORERS"]); i++){
		fColIDSum[i] = SDMan->GetCollectionID("detector2/eDep_" + IntToString(i));
	}
}

B1Run::~B1Run(){
	delete [] fColIDSum;
	delete [] fMapSum;
}

void B1Run::RecordEvent(const G4Event* evt)
{
  G4Run::RecordEvent(evt);


  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
  for (G4int i=0;i<(pp->int_map["i_NUM_OF_SCORERS"]);i++){
	  G4THitsMap<G4double>* evtMap
		= (G4THitsMap<G4double>*)(HCE->GetHC(fColIDSum[i]));
	  fMapSum[i] += *evtMap;
  }
}

void B1Run::Merge(const G4Run * aRun) {
  const B1Run * localRun = static_cast<const B1Run *>(aRun);
  for(G4int i=0;i<(pp->int_map["i_NUM_OF_SCORERS"]) ;i++){
	  fMapSum[i] += localRun->fMapSum[i];
  }
  G4Run::Merge(aRun);
}

G4double B1Run::GetTotal(const G4THitsMap<G4double> &map) const
{
  G4double tot = 0.;
  if(map.GetSize()==0) return tot;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++)
  { tot += *(itr->second); }
  return tot;
}



