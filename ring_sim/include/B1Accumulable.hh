/*
 * B1Accumulable.hh
 *
 *  Created on: Feb 20, 2018
 *      Author: adamgeva
 */

#ifndef B1ACCUMULABLE_HH_
#define B1ACCUMULABLE_HH_

#include "G4VAccumulable.hh"
#include "globals.hh"
#include "params.hh"
#include "globalFunctions.hh"
#include "GradientKey.hh"
#include <map>
#include "ParamsSingleton.hh"


class B1Accumulable : public G4VAccumulable
{
 public:
	B1Accumulable(const G4String& name);

	virtual ~B1Accumulable();

	void updateP(G4int voxel_ind);
	void updateSm_hat(G4int voxel, G4int element, G4int detector, G4double value);
	void writeGradientAndP(G4int runNum);
	virtual void Merge(const G4VAccumulable& other);
	virtual void Reset();
	void SetErrorArray(G4double* error_arr)
		{ferror_arr=error_arr;}

	private:
    // gradient array - built per thread and holds the current gradient of all replicas (detector elements) w.r.t voxels and elements.
    //G4double fSm_hat[NUM_OF_VOXELS][NUM_OF_ELEMENTS] = {};
	//std::map<GradientKey,G4double> fSm_hat;
	G4double** fSm_hat;

    //G4double fP[NUM_OF_VOXELS] = {};
	//std::map<G4int,G4double> fP;
	G4double* fP;
	G4double* ferror_arr;
	ParamsSingleton* pp = ParamsSingleton::Instance();


};



#endif /* B1ACCUMULABLE_HH_ */
