/*
 * B1EnergyDeposit.cc
 *
 *  Created on: Apr 22, 2017
 *      Author: adamgeva
 */

#include "B1EnergyDeposit.hh"
#include "globalFunctions.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"
#include "G4Run.hh"

#include <math.h>

B1EnergyDeposit::B1EnergyDeposit(G4String name, G4int type)
//per thread
: G4PSEnergyDeposit(name), fGradAccum(NULL)
{
	fscorerType =  type;
	if (type==0){ //randomly pick the 0 scorer
		// Get accumulable from the accumulable manager
		if (pp->int_map["i_CALC_GRADIENT"] == 1){
		G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
		fGradAccum = (B1Accumulable*)accumulableManager->GetAccumulable(0);
		}
	}
}

B1EnergyDeposit::~B1EnergyDeposit()
{
}

G4int B1EnergyDeposit::GetIndex(G4Step* step){
	G4StepPoint* preStepPoint = step->GetPreStepPoint();
	//entering the detector
	G4TouchableHistory* touchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4int ReplicaNum0 = touchable->GetReplicaNumber(0);
	G4int ReplicaNum1 = touchable->GetReplicaNumber(1);
	//G4int ReplicaNum2 = touchable->GetReplicaNumber(2);
	//10 is the number of bins
	return (ReplicaNum0);
//	return (ReplicaNum0 + ReplicaNum1*pp->int_map["i_NUM_OF_DET_COLS"]);

}

G4bool B1EnergyDeposit::ProcessHits(G4Step* aStep,G4TouchableHistory* touchable){
	//test: getting information from tracks
	G4Track* track = aStep->GetTrack();
//	G4double trackWeight = track->GetWeight();

	G4ThreeVector momentum_direction = track->GetMomentumDirection();
	//origin
	
	G4double x0 = -pp->double_map["f_SOURCE_TO_CENTER"]*mm;
	G4double y0 = 0;
	G4double z0 = 0;
	G4double diffx = track->GetPosition().getX()/mm;
	G4double diffy = track->GetPosition().getY()/mm;
	G4double diffz = track->GetPosition().getZ()/mm;
	
	G4double dot_prod = diffx * momentum_direction.getX()/mm + diffy * momentum_direction.getY()/mm + diffz * momentum_direction.getZ()/mm;
	G4double norm_det_angle = sqrt(pow(diffx,2) + pow(diffy,2) + pow(diffz,2));
	G4double norm_phot_angle = sqrt(pow(momentum_direction.getX(),2) + pow(momentum_direction.getY(),2) + pow(momentum_direction.getZ(),2));
	G4double cos_theta = dot_prod/(norm_det_angle*norm_phot_angle);


	//std::cout << "cos_theta " << cos_theta << std::endl;
	G4bool result = FALSE;
	G4VUserTrackInformation* info = track->GetUserInformation();
	B1TrackInformation* theInfo = (B1TrackInformation*)info;
	G4int totalNumOfInteractions = theInfo->GetNumberOfCompton() + theInfo->GetNumberOfRayl();
	G4int stepnumber = track->GetCurrentStepNumber();
	
	if(theInfo->hit_detector[fscorerType]==0)
		theInfo->hit_detector[fscorerType] = stepnumber;
	
	if (stepnumber > theInfo->hit_detector[0])
		return FALSE;
	
	
	//printing debug info:
	// G4cout << "track ID: " << track->GetTrackID() << " track parent ID: " << track->GetParentID() << " step number: " << track->GetCurrentStepNumber() << " weight " << weight << " ENERGY:  " <<   aStep->GetPreStepPoint()->GetTotalEnergy()/keV << "num of compt: " << theInfo->GetNumberOfCompton() <<  "num of rayl: " << theInfo->GetNumberOfRayl() << G4endl;

	//very rare case in biasing - in case the photon underwent an interaction but wasnt split even tho it should have, example:
	//a photon hits the detector, potoelectric absorption and the an emittion of a new photon (ID 1 for example) then this photon undergoes Rayl in the phantom and arrives at the detector with weight 1.
	//TODO: check if correct!
	//not recording the main photons - only the virtuals
	
	
	// G4cout << "SCORER: "<< fscorerType << "****ENERGY:  "  <<   aStep->GetPreStepPoint()->GetTotalEnergy()/keV << G4endl;  // aStep->GetPreStepPoint()->GetTotalEnergy()
	// G4cout << "ParentID = " << track->GetParentID()  << G4endl;
	if ((pp->int_map["i_BIASING"]==1) && (track->GetParentID()==0)) {
		result = FALSE;
	}
	else if ((pp->int_map["i_SINGLE_SCATTERING"]==1) && (totalNumOfInteractions>1)){
		result = FALSE;
	}
	else if(fscorerType==0){

		if (totalNumOfInteractions==0)
		{			
			theInfo->hit_detector[fscorerType]=stepnumber;
			G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
			pp->run_images_arr[fscorerType][GetIndex(aStep)] += edep;
			result = G4PSEnergyDeposit::ProcessHits(aStep,touchable);
		}
		else
		{
			result = FALSE;
		}
	}
	else if(fscorerType==1){

		if (totalNumOfInteractions==1)
		{			
			theInfo->hit_detector[fscorerType]=stepnumber;
			G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
			pp->run_images_arr[fscorerType][GetIndex(aStep)] += edep;
			result = G4PSEnergyDeposit::ProcessHits(aStep,touchable);
		}
		else
		{
			result = FALSE;
		}
	}
	// else if (fscorerType==1 || fscorerType==0){
		// G4cout << "****DETECTOR: " <<GetIndex(aStep) << G4endl;
		// result = recordInteraction(aStep,touchable,totalNumOfInteractions,fscorerType, theInfo);
	// }
	else if(fscorerType==2){
		theInfo->hit_detector[fscorerType]=stepnumber;
		G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
		G4int ind = GetIndex(aStep);
		pp->run_images_arr[fscorerType][ind] += edep;
		result = G4PSEnergyDeposit::ProcessHits(aStep,touchable);
	}
	else if(fscorerType==3){

		if ((theInfo->GetNumberOfCompton()<2) && (theInfo->GetNumberOfRayl()==0) && (totalNumOfInteractions==1))
		{			
			theInfo->hit_detector[fscorerType]=stepnumber;
			G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
			pp->run_images_arr[fscorerType][GetIndex(aStep)] += edep;
			result = G4PSEnergyDeposit::ProcessHits(aStep,touchable);
		}
		else
		{
			result = FALSE;
		}
	}
	else if(fscorerType==4){
		if ((theInfo->GetNumberOfRayl()<2) && (theInfo->GetNumberOfCompton()==0) && (totalNumOfInteractions==1)){
			theInfo->hit_detector[fscorerType]=stepnumber;
			G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
			pp->run_images_arr[fscorerType][GetIndex(aStep)] += edep;
			result = G4PSEnergyDeposit::ProcessHits(aStep,touchable);
		}
		else
		{
			result = FALSE;
		}
	}

	else if(fscorerType>4 && fscorerType!=16 ){
		result = recordInteraction_extra(aStep,touchable,totalNumOfInteractions,fscorerType, theInfo);
	}
	else if(fscorerType == 16){ //anti scatter grid
		result = recordInteraction_anti_scatter(aStep,touchable,cos_theta);
	}

	// this scorer is in charge of writing the path file, could be any scorer.
	if(fscorerType==0 && pp->int_map["i_CALC_GRADIENT"] == 1){
		//energy deposited in the detector - equals to the final photon energy???
		G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;

		//calc gradient:

		G4int detIndex = GetIndex(aStep);
	   	while (!theInfo->fpathLogList.empty()){
	   		segment seg = theInfo->fpathLogList.front();
	   		theInfo->fpathLogList.pop_front();
	   		if (pp->int_map["i_CALC_GRADIENT"] == 1){
	   			if (pp->int_map["i_SINGLE_SCATTERING"] == 0){
	   					updateGradTable(seg,edep,detIndex);
	   			} else if (totalNumOfInteractions < 2){
   						updateGradTable(seg,edep,detIndex);
	   			}
	   		}
	   	}
		
	}

	return result;
}

G4bool B1EnergyDeposit::recordInteraction (G4Step* aStep,G4TouchableHistory* touchable, G4int totalNumOfInteractions, G4int i, B1TrackInformation* theInfo) {
	if (totalNumOfInteractions!=i){//recording all photons that did not interact in the phantom - transmission
		return FALSE;
	} else {
		//G4cout << "recording" << i << G4endl;
		theInfo->hit_detector[fscorerType]=1;
		G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
		pp->run_images_arr[fscorerType][GetIndex(aStep)] += edep;
		return G4PSEnergyDeposit::ProcessHits(aStep,touchable);
	}
}

G4bool B1EnergyDeposit::recordInteraction_extra (G4Step* aStep,G4TouchableHistory* touchable, G4int totalNumOfInteractions, G4int i, B1TrackInformation* theInfo) {
	if (totalNumOfInteractions != (i-pp->int_map["i_NUM_OF_BASE_SCORERS"])){
		return FALSE;
	} else {
		//G4cout << "recording" << i << G4endl;
		theInfo->hit_detector[fscorerType]=1;
		G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
		pp->run_images_arr[fscorerType][GetIndex(aStep)] += edep;
		return G4PSEnergyDeposit::ProcessHits(aStep,touchable);
	}
}

G4bool B1EnergyDeposit::recordInteraction_anti_scatter (G4Step* aStep,G4TouchableHistory* touchable, G4double cos_theta) {
	G4double theta =  acos (cos_theta) * 180.0 / PI;
	G4bool survive;
	if (abs(theta) > pp->double_map["f_THETA_CUT"]){
		survive = FALSE;
	} else {
		G4double p = -(pp->double_map["f_FILL_FACTOR"]/pp->double_map["f_THETA_CUT"])
				* (abs(theta) - pp->double_map["f_THETA_CUT"]);
		G4double rndm = G4UniformRand();
		if (rndm <= p) {
			survive = TRUE;
		} else {
			survive = FALSE;
		}
	}
	if (survive){
		G4double  edep = aStep->GetPreStepPoint()->GetTotalEnergy()/keV;
		pp->run_images_arr[fscorerType][GetIndex(aStep)] += edep;
		return G4PSEnergyDeposit::ProcessHits(aStep,touchable);
	} else {
		return FALSE;
	}
}


// methods for grad calculations
G4double B1EnergyDeposit::getTotalMicXS(G4Element* el, G4double Energy){
	// currently dealing with compton and photoelectric absorption only
	G4EmCalculator emCalculator;
	G4double comptXS;
	G4double photXS;
	G4double raylXS;
	photXS = emCalculator.ComputeCrossSectionPerAtom(Energy,"gamma","phot",el,0)/cm2;
	comptXS = emCalculator.ComputeCrossSectionPerAtom(Energy,"gamma","compt",el,0)/cm2;
	raylXS = emCalculator.ComputeCrossSectionPerAtom(Energy,"gamma","Rayl",el,0)/cm2;
	return photXS + comptXS + raylXS;
}


G4double B1EnergyDeposit::getComptonMicDifferentialXS(G4Element* el, G4double E0 , G4double E1){
	G4double electron_mass_c2 = 510.99906 * keV;
	G4double classic_electr_radius = 2.818e-13 * cm;
	G4double classic_electr_radius2 = classic_electr_radius * classic_electr_radius;
	G4double eps = E1/E0;
	G4double eps_round = round( eps * 10000.0 ) / 10000.0; //rounds to 4 decimal points
	G4double deps = 1e-4; //resolution
	G4double Z = el->GetZ();
	G4double cost = 1 - electron_mass_c2/(eps_round*E0) + electron_mass_c2/E0;
	G4double cost2 = cost * cost;
	G4double sint2 = 1 - cost2;
	G4double f = 1/eps_round + eps_round;
	G4double q = 1 - (eps_round * sint2)/(1 + eps_round * eps_round);
	G4double dsigma_deps = M_PI * classic_electr_radius2 * (electron_mass_c2/E0) * f * q * Z / (2 * M_PI);
	return dsigma_deps; //* deps;
}

G4double B1EnergyDeposit::getRaylMicDifferentialXS(G4Element* el, G4double Energy, G4double angle){
	G4EmCalculator emCalculator;
	G4double cost = cos(angle);
	G4double cost2 = cost * cost;
	G4double sint = sin(angle);
	G4double fact = (3.0/8.0);
	G4double phase = fact * (1 + cost2) * sint ;
	G4double raylXS = emCalculator.ComputeCrossSectionPerAtom(Energy,"gamma","Rayl",el,0)/cm2;
	return raylXS * phase;
}

G4double B1EnergyDeposit::getComptonMacDifferentialXS(G4Material* mat, G4double E0 , G4double E1){
	const G4ElementVector* curr_element_vector = mat->GetElementVector();
	//vector of number of atoms per volume
	const G4double* curr_num_of_atoms =  mat->GetVecNbOfAtomsPerVolume();
	G4int nElements = mat->GetNumberOfElements();
	//calc dsigma macroscopic:
	G4double dsigmaMac = 0;
	for (G4int i=0 ; i<nElements ; i++) {
		dsigmaMac = dsigmaMac + curr_num_of_atoms[i]/(1/cm3) * getComptonMicDifferentialXS((*curr_element_vector)[i], E0 , E1);
	}
	return dsigmaMac;
}

G4double B1EnergyDeposit::getRaylMacDifferentialXS(G4Material* mat, G4double E0 , G4double angle){
	const G4ElementVector* curr_element_vector = mat->GetElementVector();
	//vector of number of atoms per volume
	const G4double* curr_num_of_atoms =  mat->GetVecNbOfAtomsPerVolume();
	G4int nElements = mat->GetNumberOfElements();
	//calc dsigma macroscopic:
	G4double dsigmaMac = 0;
	for (G4int i=0 ; i<nElements ; i++) {
		dsigmaMac = dsigmaMac + curr_num_of_atoms[i]/(1/cm3) * getRaylMicDifferentialXS((*curr_element_vector)[i], E0 , angle);
	}
	return dsigmaMac;
}

void B1EnergyDeposit::updateGradTable(segment seg, G4double final_energy, G4int detIndex){
	//calc da_di_dvox - gradient of the attenuation w.r.t all the elements in the voxel
	// - 1 to have the absolute voxel number in the phantom grid
	G4int curr_voxel = seg.voxel - 1;
	//we dont calculate the gradient of the air voxel
	if (curr_voxel == -1){
		return;
	}
	//update P array:
	fGradAccum->updateP(curr_voxel);
	G4double curr_energy = seg.incidentEnergy/keV;	
	G4double curr_len = seg.pathLen;
	G4double curr_density = seg.Mat->GetDensity()/(g/cm3);
	//vector of elements in current material/voxel
	const G4ElementVector* curr_element_vector = seg.Mat->GetElementVector();
	const G4double* curr_frac_vector = seg.Mat->GetFractionVector();
	//vector of number of atoms per volume
	const G4double* curr_num_of_atoms =  seg.Mat->GetVecNbOfAtomsPerVolume();
	const G4String mat_name = seg.Mat->GetName();
	char mat = mat_name[3];
	char zero = '0';
	G4int nElements = seg.Mat->GetNumberOfElements();
	//incase there is compton scatter - calc dsigma:
	
	// G4cout << "LEN = " << curr_len << " DENS = " << curr_density << " ENERGY = " << curr_energy << " INC EN = "  << seg.incidentEnergy/keV << " SCAT EN = " << seg.scatteredEnergy/keV << " FIN EN = " << final_energy << G4endl;
	
	if (seg.endingProcess == "compt"){
		//calc dsigma macroscopic:
		G4double dsigmaMac = getComptonMacDifferentialXS(seg.Mat,seg.incidentEnergy/keV,seg.scatteredEnergy/keV);
		// iterate over elements and calc attenuation factor and scatter:
		for (G4int i=0 ; i<nElements ; i++) {
				G4double n_i = curr_num_of_atoms[i]/(1/cm3);
				if (n_i == 0) continue;
				G4double frac_i = curr_frac_vector[i];
				G4double N_by_A = n_i / (curr_density * frac_i);
				G4Element* el_i =  (*curr_element_vector)[i];
				G4double dsigma_di_dvox = N_by_A * getComptonMicDifferentialXS(el_i,seg.incidentEnergy/keV,seg.scatteredEnergy/keV) * (1/dsigmaMac);
				G4double da_di_dvox = -1 * curr_len * N_by_A * getTotalMicXS(el_i,curr_energy*keV);
				//update gradient:
				fGradAccum->updateSm_hat(curr_voxel,i,detIndex,final_energy * (da_di_dvox + dsigma_di_dvox));
		}
	}
	else if (seg.endingProcess == "Rayl"){
		//calc dsigma macroscopic:
		G4double dsigmaMac = getRaylMacDifferentialXS(seg.Mat,seg.incidentEnergy/keV,seg.scatteredAngle);
		// iterate over elements and calc attenuation factor and scatter:
		for (G4int i=0 ; i<nElements ; i++) {	
				G4double n_i = curr_num_of_atoms[i]/(1/cm3);
				if (n_i == 0) continue;
				G4double frac_i = curr_frac_vector[i];
				G4double N_by_A = n_i / (curr_density * frac_i);
				G4Element* el_i =  (*curr_element_vector)[i];
				G4double dsigma_di_dvox = N_by_A * getRaylMicDifferentialXS(el_i,seg.incidentEnergy/keV,seg.scatteredAngle) * (1/dsigmaMac);
				G4double da_di_dvox = -1 * curr_len * N_by_A * getTotalMicXS(el_i,curr_energy*keV);
				//update gradient:
				fGradAccum->updateSm_hat(curr_voxel,i,detIndex,final_energy * (da_di_dvox + dsigma_di_dvox));
		}
	}
	else { //no compton or rayl at the end
		for (G4int i=0 ; i<nElements ; i++) {
				G4double n_i = curr_num_of_atoms[i]/(1/cm3);
				G4double frac_i = curr_frac_vector[i];
				if (n_i == 0) continue;

				G4double N_by_A = n_i / (curr_density * frac_i);
				G4Element* el_i =  (*curr_element_vector)[i];
				std::string el_name = el_i->GetName();
				//G4double micmic = getTotalMicXS(el_i,curr_energy*keV);
				G4double da_di_dvox = -1 * curr_len * N_by_A * getTotalMicXS(el_i,curr_energy*keV);
				//update gradient:
				fGradAccum->updateSm_hat(curr_voxel,i,detIndex,final_energy * (da_di_dvox));

		}
	}
}
