
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B1StackingAction.hh"


#include "G4RunManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1StackingAction::B1StackingAction()
 : G4UserStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1StackingAction::~B1StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
B1StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{

  //keep primary particle
	
   // G4cout << "PARENTID = " << aTrack->GetParentID() << G4endl;
  if (aTrack->GetParentID() == 0) { return fUrgent; }

  //stack or delete secondary electrons
  G4ClassificationOfNewTrack status = fUrgent;
  if (aTrack->GetDefinition()->GetParticleName() == "e-") {
	 // G4cout << "killed electron" << G4endl;
	status = fKill;
  }


  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
