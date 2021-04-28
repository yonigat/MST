#include "B1PhantomParameterisationColour.hh"

#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "params.hh"

B1PhantomParameterisationColour::B1PhantomParameterisationColour()
: G4PhantomParameterisation()
{

    //ReadColourData();
	if (pp->int_map["i_CALC_GRADIENT"] == 1){
		SetSkipEqualMaterials(false);
	}
	else {
		SetSkipEqualMaterials(true);
	}
}

B1PhantomParameterisationColour::~B1PhantomParameterisationColour()
{
}

void B1PhantomParameterisationColour::ReadColourData()
{
    //----- Add a G4VisAttributes for materials not defined in file;
    
}

G4Material* B1PhantomParameterisationColour::
ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol, const G4VTouchable *)
{
    G4Material* mate = G4PhantomParameterisation::ComputeMaterial( copyNo, physVol, 0 );
    
    return mate;
}
