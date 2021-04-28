#include "params.hh"
#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "globalFunctions.hh"
#include <iostream>
#include <cmath>
#include <time.h>
#include <math.h>
#include <map>
#include <bits/random.h>
#include <string>
#include <sstream>
#include <fstream>



B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
	this->fParticleGun = new G4ParticleGun[pp->int_map["i_NUM_OF_SOURCES"]]; // pointer a to G4 gun class array, containing all the sources around the ring
	this->code = new G4int*[pp->int_map["i_NUM_SHOTS"]];
	this->intens_code = new G4int*[pp->int_map["i_NUM_SHOTS"]];
	for (int i = 0; i<pp->int_map["i_NUM_SHOTS"]; i++ ) {
     		code[i] = new G4int[pp->int_map["i_NUM_PROJECTION_IN_SHOT"]];
     		intens_code[i] = new G4int[pp->int_map["i_NUM_PROJECTION_IN_SHOT"]];
    }
	sum_dose = 0;
	readSpectrum();
	readCode();
	G4int n_particle = 1;
	// default particle kinematic
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName =
			"gamma");
	int iDims = 3;
	G4double og = 1*mm;
	G4double loc_x, loc_y, loc_z;
	G4double curr_loc;
	for (int i=0; i<pp->int_map["i_NUM_OF_SOURCES"]; i++){
		for (int j=0; j<iDims; j++) {
			G4int loc_copy = j + i*iDims;
			curr_loc = pp->src_arr[loc_copy];
			if (j==0) loc_x = curr_loc;
			else if (j==1) loc_y = curr_loc;
			else if (j==2) loc_z = curr_loc;
		}
		G4ThreeVector loc = G4ThreeVector(loc_x, loc_y, loc_z)*og;
		fParticleGun[i].SetNumberOfParticles(n_particle);
		fParticleGun[i].SetParticleDefinition(particle);
		fParticleGun[i].SetParticleEnergy(pp->double_map["f_PARTICLE_ENERGY"] * keV);
		fParticleGun[i].SetParticlePosition(loc);
	}
}

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction() {
	delete [] fParticleGun;
	for (int i = 0; i<pp->int_map["i_NUM_SHOTS"]; i++ ) {
     		delete [] code[i];
     		delete [] intens_code[i];
    }
	delete [] code;
	delete [] intens_code;
	std::cout << "****** Overall dosage:" << sum_dose << std::endl;
	pp->num_phot = 0;
}

void B1PrimaryGeneratorAction::readCode() {

	G4int curr_source = 0;
	G4int curr_intens = 0;
	// Setup the weights (in this case linearly weighted)
	for (int i = 0; i < pp->int_map["i_NUM_SHOTS"]; i++) {
		for (int j = 0 ; j < pp->int_map["i_NUM_PROJECTION_IN_SHOT"]; j++){
			G4int code_copy = j + i*pp->int_map["i_NUM_PROJECTION_IN_SHOT"];
			curr_source = pp->projection_code_arr[code_copy];
			curr_intens = pp->intens_code_arr[code_copy];
			code[i][j] = curr_source;
			intens_code[i][j] = curr_intens;
		}
	}
}

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4int* runID = new G4int[pp->int_map["i_NUM_PROJECTION_IN_SHOT"]];
	G4int* intens_runID = new G4int[pp->int_map["i_NUM_PROJECTION_IN_SHOT"]];
		fAllSources = 1;
		pp->num_phot++;
		num_phot = pp->num_phot;
		if (fAllSources == 0) {
			//runID = frunIDRand;
		} else if (fAllSources == 1){
			G4int x = G4RunManager::GetRunManager()->GetNonConstCurrentRun()->GetRunID();
			for (int k = 0 ; k < pp->int_map["i_NUM_PROJECTION_IN_SHOT"]; k++){
			runID[k] = code[x][k];
			intens_runID[k] = intens_code[x][k];
			}
		} else {
		//	runID =
		//			G4RunManager::GetRunManager()->GetNonConstCurrentRun()->GetRunID();
		}

		G4double u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z;
		G4double src2det_dist = pp->double_map["f_SOURCE_TO_CENTER"]*mm + pp->double_map["f_CENTER_TO_DET"]*mm;

		//		G4double delta_dist =
		for (int l = 0; l < pp->int_map["i_NUM_PROJECTION_IN_SHOT"]; l++){
			G4double curr_orient;
			G4int iDims = 3;
			G4int curr_source = runID[l];
			G4int curr_intens = intens_runID[l];
			// std::cout << "curr_intens:     " << curr_intens << std::endl;
			if (num_phot > curr_intens) {
				// std::cout << num_phot << ">" <<  curr_intens << std::endl;
				continue;
			}
			for (int iy = 0; iy <iDims; iy++) {
				for (int ix = 0 ; ix <iDims; ix++){
					G4int code_copy = ix + iy*iDims + curr_source*iDims*iDims;
					curr_orient = pp->src_orient_arr[code_copy];
					if (iy == 0 && ix == 0) u_x = curr_orient;
					else if (iy == 0 && ix == 1) v_x = curr_orient;
					else if (iy == 0 && ix == 2) w_x = curr_orient;
					else if (iy == 1 && ix == 0) u_y = curr_orient;
					else if (iy == 1 && ix == 1) v_y = curr_orient;
					else if (iy == 1 && ix == 2) w_y = curr_orient;
					else if (iy == 2 && ix == 0) u_z = curr_orient;
					else if (iy == 2 && ix == 1) v_z = curr_orient;
					else if (iy == 2 && ix == 2) w_z = curr_orient;
				}
			}

			G4ThreeVector theta_vec = G4ThreeVector(u_x, u_y, u_z);
			G4ThreeVector delta_vec = G4ThreeVector(v_x, v_y, v_z);
			G4ThreeVector prop_vec = G4ThreeVector(w_x, w_y, w_z);
			G4RotationMatrix rotm  = G4RotationMatrix(prop_vec, theta_vec, delta_vec);

			G4double MinTheta = pp->double_map["f_MIN_THETA"];
			G4double MaxTheta = pp->double_map["f_MAX_THETA"];
			G4double MinPhi = pp->double_map["f_MIN_PHI"];
			G4double MaxPhi = pp->double_map["f_MAX_PHI"];

			G4ThreeVector total_vec;
			G4double phi;

				//phi

			if (pp->int_map["i_CONE_BEAM"] == 1){
				G4double rndm1 = G4UniformRand();
				phi = MinPhi + rndm1 * (MaxPhi - MinPhi);
				//cos,sin theta - used for cone beam
				G4double rndm = G4UniformRand();
				G4double costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta) - std::cos(MaxTheta));
				G4double sintheta = std::sqrt(1. - costheta*costheta);
				//cos,sin theta - used for fan beam
				//cos,sin phi
				G4double cosphi = std::cos(phi);
				G4double sinphi = std::sin(phi);
				//coordinates
				G4double px = costheta;
				G4double py = sintheta * cosphi;
				G4double pz = sintheta * sinphi;
				total_vec = G4ThreeVector(px, py, pz);
				// total_vec = G4ThreeVector(1.0, 0, 0);
			} else {

				//cos,sin theta - used for cone beam
				G4double rndm = G4UniformRand();
				G4double theta =(MaxTheta)*(2*rndm-1);
				// G4double theta = std::acos((std::cos(MaxTheta+M_PI/2)*(2*rndm-1)))-M_PI/2;
				G4double costheta = std::cos(theta);
				G4double sintheta = std::sin(theta);
				//cos,sin phi
				G4double cosphi = std::cos(phi);
				G4double sinphi = std::sin(phi);
				//coordinates
				G4double px = costheta;
				G4double py = sintheta * cosphi;
				total_vec = G4ThreeVector(px, py, 0.000001);

			}


			G4double energyP;
			G4int energyBin;

			//if switching sources is applied

			if (pp->int_map["i_ALT_SOURCES"] == 1) {
				// sourceAngleDiff is the angle between every source
				total_vec = rotm * total_vec;
				fParticleGun[runID[l]].SetParticleMomentumDirection(total_vec);
				energyBin = getEnergyInd();//+1;
				energyP = energyBin * keV;
				sum_dose += energyBin;
				G4ThreeVector pos = fParticleGun[runID[l]].GetParticlePosition();
				fParticleGun[runID[l]].SetParticleEnergy(energyP);
				// fParticleGun[0].SetParticleEnergy(50*keV);
				fParticleGun[runID[l]].GeneratePrimaryVertex(anEvent);
			} else {
				fParticleGun[0].SetParticleMomentumDirection(total_vec);
				energyBin = getEnergyInd();
				sum_dose += energyBin;
				energyP = energyBin * keV;
				fParticleGun[0].SetParticleEnergy(energyP);
//				fParticleGun[0].SetParticleEnergy(120*keV);
				fParticleGun[0].GeneratePrimaryVertex(anEvent);
			}
		}

	delete [] runID;
	delete [] intens_runID;

}

void B1PrimaryGeneratorAction::readSpectrum() {
//	G4double spectrum_dist[pp->int_map["i_NUM_OF_SPECTRUM_BINS]]; //distribution of energy for every source
	G4double num_of_photons = 0;

	// Setup the weights (in this case linearly weighted)
	for (int i = 0; i < pp->int_map["i_NUM_OF_SPECTRUM_BINS"]; i++) {

		num_of_photons = pp->spectrum_arr[i];
		
		// std::cout << "i= " << i << " spect = " <<  num_of_photons << std::endl;
//			spectrum_dist[i] = ceil(num_of_photons * NUM_OF_PHOTONS);
		fweights.push_back((num_of_photons));
	}
//	G4int sum = 0;
//	for (int k = 0; k<NUM_OF_SPECTRUM_BINS; k ++){
//		std::cout << "i= " << k << " spect = " <<  fspectrum_dist[k][0] << std::endl;
//		sum = sum + fspectrum_dist[k][0];
//	}
}

G4int B1PrimaryGeneratorAction::getEnergyInd(){
	//todo - maybe we can declare discrete_distribution globally to avoid initialization
	std::discrete_distribution<G4int> distribution(fweights.begin(),fweights.end());
	std::random_device rd;
	fgenerator.seed(rd());
	G4int bin = distribution(fgenerator);
	return bin;

//	std::discrete_distribution<G4int> distribution(fweights.begin(), fweights.end());
//	G4int bin = distribution(fgenerator);
//	return bin;
}




