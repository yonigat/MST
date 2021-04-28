#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4Box.hh"
#include "params.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
class G4Material;

#include <vector>
#include <set>
#include <map>
#include "ParamsSingleton.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;

struct materialIDs {
  G4int mat1ID;
  G4int mat2ID;
};

struct matInfo
{
  G4double fSumdens;
  G4int fNvoxels;
  G4int fId;
};

// Detector construction class to define materials and geometry.
class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    virtual void ConstructSDandField();
    
    void setContainerRotation(G4double delta);

   	std::vector<G4Material*> fMaterials;
   	size_t* fMateIDs; // index of material of each voxel - this array is in the size of the number of voxels


  protected:
    // create the original materials
    void ReadPhantomData();
    void ReadPhantomDataCT();

    void InitialisationOfMaterialsCT();
    void InitialisationOfMaterialsCT_basic();
    void InitialisationOfMaterials();

    void ReadVoxelDensities();

    void ConstructPhantomContainer();
    virtual void ConstructPhantom() = 0;
	// construct the phantom volumes.
	//  This method should be implemented for each of the derived classes

    //*************************************XCAT
    void InitialisationOfMaterialsMap_XCAT();
    void InitialisationOfMaterials_XCAT();
    void ReadPhantomData_XCAT();
    void ReadPhantomDataFile_XCAT(const G4String& fname, G4int sliceNumber);

    G4Material* BuildMaterialWithChangingDensity(const G4Material* origMate, float density, G4String newMateName );



  protected:
    G4LogicalVolume* fvoxel_logic;
    G4Box* fContainer_solid;
   	G4LogicalVolume* fContainer_logic;
   	G4VPhysicalVolume* fContainer_phys;

   	//std::vector<G4Material*> fMaterials;
   	//size_t* fMateIDs; // index of material of each voxel - this array is in the size of the number of voxels

   	G4int fNVoxelX, fNVoxelY, fNVoxelZ;
	G4double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;
	G4double fMinX,fMinY,fMinZ; // minimum extension of voxels (position of wall)
	G4double fMaxX,fMaxY,fMaxZ; // maximum extension of voxels (position of wall)
	G4RotationMatrix phantom_rotm;
	G4ThreeVector phantom_loc;
    std::map<G4int,G4Material*> thePhantomMaterialsOriginal;
     // map numberOfMaterial to G4Material. They are the list of materials as built from .geom file
	std::map<G4int,G4Material*> fBaseMaterials;
	// map numberOfMaterial to G4Material. They are the list of materials as built from .geom file

	//***************************************XCAT
	std::map<G4int,materialIDs> intensityToMateID;
	// maps the intensity value in every voxel to its material ID


  private:
	G4LogicalVolume* worldLV;
	//logical detectors
    G4LogicalVolume* detectorLV;
    std::vector<G4VisAttributes*> fVisAttributes;

   	ParamsSingleton* pp = ParamsSingleton::Instance();


};


#endif

