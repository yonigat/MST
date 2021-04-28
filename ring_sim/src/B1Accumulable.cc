/*
 * B1Accumulable.cc
 *
 *  Created on: Feb 20, 2018
 *      Author: adamgeva
 */



#include "B1Accumulable.hh"
#include "G4VAccumulable.hh"
#include "params.hh"

#include "globalFunctions.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include "G4ios.hh"
#include "globals.hh"



B1Accumulable::B1Accumulable(const G4String& name):
G4VAccumulable(name),ferror_arr(0)
{
	fSm_hat = new G4double*[pp->int_map["i_NUM_OF_VOXELS"]];

	for(int i = 0; i < pp->int_map["i_NUM_OF_VOXELS"]; ++i)
		fSm_hat[i] = new G4double[pp->int_map["i_NUM_OF_ELEMENTS"]];

	fP= new G4double[pp->int_map["i_NUM_OF_VOXELS"]];

}

B1Accumulable:: ~B1Accumulable() {
	for(int i = 0; i < pp->int_map["i_NUM_OF_VOXELS"]; ++i)
			delete [] fSm_hat[i];
	delete [] fSm_hat;
	delete [] fP;

}

void B1Accumulable::updateP(G4int voxel_ind){
	//update P array:
	fP[voxel_ind]++;
}

void B1Accumulable::updateSm_hat(G4int voxel, G4int element, G4int detector, G4double value){
	G4double error = ferror_arr[detector];
	// G4cout << "***DETECTOR = " << detector << " ERROR = " << ferror_arr[detector] << G4endl;
	//update gradient
	fSm_hat[voxel][element] = fSm_hat[voxel][element] + value * error;
	//GradientKey curr_key = GradientKey(voxel,element);
	//fSm_hat[curr_key] = fSm_hat[curr_key] + value * error;
}

void B1Accumulable::Merge(const G4VAccumulable& other){
	const B1Accumulable& otherAccumulable
	= static_cast<const B1Accumulable&>(other);

	//iterate and accumulate
	for (int i=0; i<pp->int_map["i_NUM_OF_VOXELS"]; i++){
		fP[i] += otherAccumulable.fP[i];
		for (int j=0; j<pp->int_map["i_NUM_OF_ELEMENTS"]; j++){
			fSm_hat[i][j] += otherAccumulable.fSm_hat[i][j];
		}
	}

	//accumulate fP:
/*
	for (std::map<G4int,G4double>::const_iterator it=otherAccumulable.fP.begin(); it!=otherAccumulable.fP.end(); ++it)
		fP[it->first] += it->second;

	//accumulate fSm_hat:
	for (std::map<GradientKey,G4double>::const_iterator it=otherAccumulable.fSm_hat.begin(); it!=otherAccumulable.fSm_hat.end(); ++it){
		fSm_hat[it->first] += it->second;
	}
*/
}


void B1Accumulable::Reset(){
	//iterate and reset arrays
	for (int i=0; i<pp->int_map["i_NUM_OF_VOXELS"]; i++){
		fP[i] = 0;
		for (int j=0; j<pp->int_map["i_NUM_OF_ELEMENTS"]; j++){
			fSm_hat[i][j] = 0;
		}
	}
	/*
	fSm_hat.clear();
	fP.clear();
	*/
	ferror_arr=NULL; // run action is in charge of deleting
}


void B1Accumulable::writeGradientAndP(G4int runNum){
	//writing gradient table to file
	if (pp->int_map["i_WRITE_TO_FILES"] == 0){
		pp->write_grad_run(runNum, fSm_hat);
		pp->write_P_run(runNum, fP);
	}else{
		std::string fileName =  std::string(pp->string_map["s_GRADIENT_DIR"])  + "grad/" + IntToString(runNum) + "run_gradient.csv";
		std::ofstream outputPathsFile_grad;
		outputPathsFile_grad.open(fileName.c_str());

		for (G4int voxel = 0; voxel < pp->int_map["i_NUM_OF_VOXELS"]; voxel ++){
			for (G4int element = 0; element < pp->int_map["i_NUM_OF_ELEMENTS"]; element ++){
				if (element < pp->int_map["i_NUM_OF_ELEMENTS"] - 1){
					G4double val = fSm_hat[voxel][element];
					outputPathsFile_grad << fSm_hat[voxel][element] << ',';
				}else
					outputPathsFile_grad << fSm_hat[voxel][element];
			}
			outputPathsFile_grad << '\n';
		}
	/*
		for (std::map<GradientKey,G4double>::iterator it=fSm_hat.begin(); it!=fSm_hat.end(); ++it){
			//writes (key,value)
			outputPathsFile_grad << it->second;
			++it;
			outputPathsFile_grad << ',' << it->second << '\n';
		}
	*/
		outputPathsFile_grad.close();

		//write P array:
		fileName =  std::string(pp->string_map["s_GRADIENT_DIR"]) + "P/" + IntToString(runNum) + "run_P.csv";
		std::ofstream outputPathsFile_P;
		outputPathsFile_P.open(fileName.c_str());

		for (G4int voxel = 0; voxel < pp->int_map["i_NUM_OF_VOXELS"]; voxel ++){
			outputPathsFile_P <<  fP[voxel] << '\n';
		}

	/*
		for (std::map<G4int,G4double>::iterator it=fP.begin(); it!=fP.end(); ++it){
			//writes (key,value)
			outputPathsFile_P <<  it->first << ',' << it->second << '\n';
		}
	*/
		outputPathsFile_P.close();
	}
}




