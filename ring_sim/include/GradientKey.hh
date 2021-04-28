/*
 * GradientKey.hh
 *
 *  Created on: Feb 28, 2018
 *      Author: adamgeva
 */

#ifndef GRADIENTKEY_HH_
#define GRADIENTKEY_HH_
class GradientKey {

public:
	GradientKey(G4int v,G4int el)
	:voxel(v), element(el){}
	G4int Getvoxel() const {return voxel;}
	G4int Getelement() const {return element;}

	bool operator< (const GradientKey& GradientKeyObj) const
	{
		if(GradientKeyObj.voxel > this->voxel) return true;
		else if (GradientKeyObj.voxel < this->voxel) return false;
		else { //equal element
			if(GradientKeyObj.element > this->element) return true;
			else if (GradientKeyObj.element < this->element) return false;
			else { //equal voxel
				return false;
			}
		}

	}

private:
	G4int voxel;
	G4int element;
};


#endif /* GRADIENTKEY_HH_ */
