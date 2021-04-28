/*
 * globalFunctions.hh
 *
 *  Created on: Apr 17, 2017
 *      Author: adamgeva
 */

#ifndef GLOBALFUNCTIONS_HH_
#define GLOBALFUNCTIONS_HH_

#define PI 3.14159265

#include "G4Types.hh"
#include <string>
#include <sstream>
#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"

//todo: remove from here to a macro!!!
std::string IntToString (int a);

G4bool RR (G4double p);

G4bool outOfRing (G4ThreeVector position, G4ThreeVector momentumDirection, G4double Zup , G4double Zdown, G4double ringRadius);

G4double angleBetweenVecs(G4ThreeVector vecA, G4ThreeVector vecB);


#endif /* GLOBALFUNCTIONS_HH_ */
