/*
 * globalFunctions.cc

 *
 *  Created on: Apr 17, 2017
 *      Author: adamgeva
 */

#include "globalFunctions.hh"
#include <stdlib.h>
#include "params.hh"
#include <math.h>




//todo: remove from here to a macro!!!
std::string IntToString (int a)
{
    std::ostringstream temp;
    temp<<a;
    return temp.str();
}
//russian roulette
//p is probability to survive. return true if particle survives and false otherwise
G4bool RR (G4double p)
{
	G4double sample = rand() / (RAND_MAX + 1.);
	if (sample < p) return true;
	return false;
}

//Check in or out of ring - true if out
G4bool outOfRing (G4ThreeVector position, G4ThreeVector momentumDirection, G4double Zup , G4double Zdown, G4double ringRadius) {
	G4double P1x = momentumDirection[0];
	G4double P1y = momentumDirection[1];
	G4double P1z = momentumDirection[2];
	G4double P0x = position[0];
	G4double P0y = position[1];
	G4double P0z = position[2];
	//normalize P0 to unit vector in x,y - p0 (small letters)
	G4double P0Norm = sqrt(pow(P0x,2)+pow(P0y,2));
	G4double p0x = P0x/P0Norm;
	G4double p0y = P0y/P0Norm;
	//normalize P1 to unit vector in x,y - p1 (small letters)
	G4double P1Norm = sqrt(pow(P1x,2)+pow(P1y,2));
	G4double p1x = P1x/P1Norm;
	G4double p1y = P1y/P1Norm;
	//calc angles
	G4double alpha = acos(-(p0x*p1x + p0y*p1y)); //radians
	G4double gamma = asin((P0Norm/ringRadius)*sin(alpha)); //radians
	//calc ll
	G4double ll = sin(alpha+gamma)*(P0Norm/sin(gamma));

	//find intersection on the ring
	G4double xRing = P0x + ll*p1x;
	G4double yRing = P0y + ll*p1y;
	//up ring vector - from current photon position to up ring
	G4double zRing = Zup;
	G4double vecx = xRing - P0x;
	G4double vecy = yRing - P0y;
	G4double vecz = zRing - P0z;
	G4double upAngle = atan(fabs(vecz)/(sqrt(pow(vecx,2)+pow(vecy,2)))) * 180 / PI;
	if (vecz<0) {upAngle=upAngle*-(-1);}
	//down ring vector
	zRing = Zdown;
	vecz = zRing - P0z;
	G4double downAngle = atan(fabs(vecz)/(sqrt(pow(vecx,2)+pow(vecy,2)))) * 180 / PI;
	if (vecz<0) {downAngle=downAngle*(-1);}
	//calc current angle
	G4double currentAngle = atan(fabs(P1z)/(sqrt(pow(P1x,2)+pow(P1y,2)))) * (180 / PI);
	if (P1z<0) {currentAngle=currentAngle*(-1);}
	//check if out bound or not
	if (currentAngle>upAngle || currentAngle<downAngle) {
		return true; //out
	}
	return false;

}

G4double angleBetweenVecs(G4ThreeVector vecA, G4ThreeVector vecB) {

	G4double P1x = vecA[0];
	G4double P1y = vecA[1];
	G4double P1z = vecA[2];
	G4double P0x = vecB[0];
	G4double P0y = vecB[1];
	G4double P0z = vecB[2];
	//normalize P0 to unit vector in x,y - p0 (small letters)
	G4double P0Norm = sqrt(pow(P0x,2)+pow(P0y,2)+pow(P0z,2));
	G4double p0x = P0x/P0Norm;
	G4double p0y = P0y/P0Norm;
	G4double p0z = P0z/P0Norm;
	//normalize P1 to unit vector in x,y - p1 (small letters)
	G4double P1Norm = sqrt(pow(P1x,2)+pow(P1y,2)+pow(P1z,2));
	G4double p1x = P1x/P1Norm;
	G4double p1y = P1y/P1Norm;
	G4double p1z = P1z/P1Norm;
	//calc angles

//	G4double alpha = acos((P0x*P1x + P0y*P1y + P0z*P1z)/(P0Norm*P1Norm)); //radians

	G4double alpha = acos(-(p0x*p1x + p0y*p1y + p0z*p1z)); //radians
	return (M_PI - alpha);
//	return (alpha);
}
//G4ThreeVector Mom = actualParticleChange->GetProposedMomentumDirection();
//G4cout << "Mom = " << Mom << G4endl;
//G4cout << "MomX = " << Mom[0] << " MomY = " << Mom[1] << " MomZ= " << Mom[2] << G4endl;
//G4cout << "Mom RMS = " << sqrt(pow(Mom[0],2) + pow(Mom[1],2) + pow(Mom[2],2)) << G4endl;
