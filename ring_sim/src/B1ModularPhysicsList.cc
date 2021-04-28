/*
 * B1ModularPhysicsList.cc

 *
 *  Created on: Sep 26, 2017
 *      Author: adamgeva
 */

#include "G4Types.hh"
#include "B1ModularPhysicsList.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "QBBC.hh"
#include "FTFP_BERT.hh"
#include "G4PhysListFactory.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4RegionStore.hh"
#include <iostream>



B1ModularPhysicsList::B1ModularPhysicsList(G4String name) :  G4VModularPhysicsList(){

	defaultCutValue  = 2.*mm;
    fCutForGamma     = defaultCutValue;


	// Physics list
	G4PhysListFactory factory;
	G4VModularPhysicsList* phys = NULL;
	phys = factory.GetReferencePhysList("FTFP_BERT_LIV");
	//phys = factory.GetReferencePhysList("FTFP_BERT");
	//FTFP_BERT_PEN
	//G4VModularPhysicsList* physicsList = new QBBC;

    if (!phys) G4Exception("B1ModularPhysicsList::B1ModularPhysicsList","InvalidSetup",
                              FatalException,"PhysicsList does not exist");

    for (G4int i = 0; ; ++i) {
       G4VPhysicsConstructor* elem =
                  const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
       if (elem == NULL) break;
       G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
       RegisterPhysics(elem);
    }
}

B1ModularPhysicsList::~B1ModularPhysicsList(){

}


void B1ModularPhysicsList::SetCuts()
{

	  SetCutsWithDefault();
	  // Production thresholds for detector regions
	  G4String regName = "ContainerRegion";
	  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regName);
	  G4ProductionCuts* cuts = new G4ProductionCuts;
	  //cuts->SetProductionCut(1*m); // same cuts  for gamma e- and e+
	  //cuts->SetProductionCut(1000*km);
	  cuts->SetProductionCut(1*cm,G4ProductionCuts::GetIndex("gamma"));
	  cuts->SetProductionCut(5000*m,G4ProductionCuts::GetIndex("e-"));
	  cuts->SetProductionCut(5000*m,G4ProductionCuts::GetIndex("e+"));
	  region->SetProductionCuts(cuts);

//	  G4ProductionCuts* cuts2 = region->GetProductionCuts();
//	  G4int index = cuts2->GetIndex("gamma");
//	  G4double pCutGamma = cuts2->GetProductionCut(index);
//	  std::cout << "gamma cut = " << pCutGamma << std::endl;
}


