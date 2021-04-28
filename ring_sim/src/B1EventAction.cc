//
#include "params.hh"
#include "B1EventAction.hh"
#include "detectorHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


B1EventAction::B1EventAction()
: G4UserEventAction(),
  fdetectorHCID(-1)
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(0);
}


B1EventAction::~B1EventAction()
{}


void B1EventAction::BeginOfEventAction(const G4Event*)
{
	//detector 1 is created only for recording histograms
	if (pp->int_map["i_RECORD_HIST"] == 1){
		if (fdetectorHCID==-1) {
		  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
		  //get the hits collection of our detector in this event
		  fdetectorHCID = sdManager->GetCollectionID("detector1/detectorColl");
		}
	}

}


void B1EventAction::EndOfEventAction(const G4Event* event)
{
	if(pp->int_map["i_RECORD_HIST"] == 1){
		// Step 1: Get the hits collection of this event
		G4HCofThisEvent* hce = event->GetHCofThisEvent();
		if (!hce)
		{
			G4ExceptionDescription msg;
			msg << "No hits collection of this event found.\n";
			G4Exception("B1EventAction::EndOfEventAction()",
						"Code001", JustWarning, msg);
			return;
		}

		// Step 2: Using the memorised IDs get the collections
		// corresponding to the detector
		// Get hits collections
		detectorHitsCollection* hHC1 = static_cast<detectorHitsCollection*>(hce->GetHC(fdetectorHCID));

		if ( (!hHC1) )
		{
			G4ExceptionDescription msg;
			msg << "Some of hits collections of this event not found.\n";
			G4Exception("B1EventAction::EndOfEventAction()",
						"Code001", JustWarning, msg);
			return;
		}
		/*
	   // Fill histograms
	   // Get analysis manager
	   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	   G4int n_hit = hHC1->entries();
	   //G4cout << "number of entries (hits) in this event: " << n_hit << G4endl;
	   for (G4int i=0;i<n_hit;i++)
	   {
		  detectorHit* hit = (*hHC1)[i];
		  G4ThreeVector localPos = hit->GetLocalPos();
		  G4double totalEnergy = hit->GetTotalEnergy();
		  //G4cout << "filling histogram with " << localPos.x()/cm << " and "<< localPos.y()/cm << " and Energy:" << totalEnergy/keV << G4endl;
		  analysisManager->FillH3(hit->GetID(), localPos.x()/cm, localPos.y()/cm,totalEnergy/keV);

	   }

		// Print diagnostics

		G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
		if ( printModulo==0 || event->GetEventID() % printModulo != 0) return;

		G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary(0);
		G4cout << G4endl
			   << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
			   << primary->GetG4code()->GetParticleName()
			   << " " << primary->GetMomentum() << G4endl;
		*/
		// Loop on the collection and dump on screen hits
		/*
		G4cout << "detector 1 has " << n_hit << " hits." << G4endl;
		for (G4int i=0;i<n_hit;i++)
		{
			detectorHit* hit = (*hHC1)[i];
			hit->Print();
		}
		*/
	}

}

