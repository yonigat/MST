#ifndef B1PhantomParameterisationColour_HH
#define B1PhantomParameterisationColour_HH

#include <map>

#include "G4PhantomParameterisation.hh"
#include "ParamsSingleton.hh"

class G4VisAttributes;

// *********************************************************************
/// Class inherited from G4PhantomParameterisation to provide different 
//  colour for each material
// *********************************************************************

class B1PhantomParameterisationColour : public G4PhantomParameterisation
{
public:  // with description
  
  B1PhantomParameterisationColour();
  ~B1PhantomParameterisationColour();
  
  virtual G4Material* ComputeMaterial(const G4int repNo, 
                                      G4VPhysicalVolume *currentVol,
                                      const G4VTouchable *parentTouch=0);
  
private:
  void ReadColourData();

private:
  std::map<G4String,G4VisAttributes*> fColours;
  ParamsSingleton* pp = ParamsSingleton::Instance();

};


#endif
