#include "B1DetectorConstruction.hh"
#include "myDetectorSD.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "globalFunctions.hh"
#include "B1EnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "B1BOptrMultiParticleChangeCrossSection.hh"
#include "B1BOptrFS.hh"
#include "B1BOptrFD.hh" //Forced Detection


#include "G4BOptrForceCollision.hh"
#include "G4EmCalculator.hh"


#include "G4UserLimits.hh"

#include <math.h>
#include <iostream>



B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fvoxel_logic(0),
  fContainer_solid(0),
  fContainer_logic(0),
  fContainer_phys(0),
  fMateIDs(0),
  fNVoxelX(0),
  fNVoxelY(0),
  fNVoxelZ(0),
  fVoxelHalfDimX(0),
  fVoxelHalfDimY(0),
  fVoxelHalfDimZ(0),
  fMinX(0),
  fMinY(0),
  fMinZ(0), // minimum extension of voxels (position of wall)
  fMaxX(0),
  fMaxY(0),
  fMaxZ(0),
  worldLV(),
  detectorLV(0),
  fVisAttributes()
{ }

B1DetectorConstruction::~B1DetectorConstruction()
{
	for (G4int i=0; i<G4int(fVisAttributes.size()); ++i)
	    {
	      delete fVisAttributes[i];
	    }
}


G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
	//set tolerance
	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*pp->double_map["f_WORLD_XY"]*cm);
	std::cout << "GetSurfaceTolerance() = " << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm << std::endl;
	std::cout << "GetAngularTolerance() = " << G4GeometryTolerance::GetInstance()->GetAngularTolerance()/mm << std::endl;
	std::cout << "GetRadialTolerance() = " << G4GeometryTolerance::GetInstance()->GetRadialTolerance()/mm << std::endl;


	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = true;

	// World
	// sizes are half size
	 // vacuum
	G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	G4double world_sizeXY = pp->double_map["f_WORLD_XY"]*cm;
	G4double world_sizeZ  = pp->double_map["f_WORLD_Z"]*cm;
	G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
	
	G4Box* worldS = new G4Box("World",world_sizeXY, world_sizeXY, world_sizeZ);
	worldLV = new G4LogicalVolume(worldS, vacuum, "World");
	G4VPhysicalVolume* worldPHS = new G4PVPlacement(0, G4ThreeVector(),worldLV,"World",0,false,0,checkOverlaps);

	// detector - specs
	G4Material* detectorMat = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	G4double detector_sizeX = pp->double_map["f_DETECTOR_X"]*mm;
	G4double detector_sizeY = pp->double_map["f_DETECTOR_Y"]*mm;
	G4double detector_sizeZ = pp->double_map["f_DETECTOR_Z"]*mm;
//	G4double dist = pp->double_map["f_CENTER_TO_DET"]*mm + pp->double_map["f_DETECTOR_X"]*mm; // detector radius + half

	// detector1
	G4Box* detectorS = new G4Box("detector",detector_sizeX, detector_sizeY, detector_sizeZ);
	detectorLV = new G4LogicalVolume(detectorS, detectorMat,"detectorLV");

	G4int num_detectors = pp->int_map["i_NUM_OF_DET_ROWS"] * pp->int_map["i_NUM_OF_DET_COLS"];
	G4int iDims = 3;
	G4double loc_x, loc_y, loc_z, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z;

	//---- Read the sources location and orientations - this is going to be ugly!
	G4double og = 1*mm;
	for( G4int iz = 0; iz < num_detectors; iz++ ) {
		for( G4int iy = 0; iy < iDims; iy++ ) {
			G4int loc_copy = iy + iz*iDims;
			G4double curr_loc = pp->det_arr[loc_copy];
			if (iy == 0) loc_x = curr_loc;
			else if (iy == 1) loc_y = curr_loc;
			else if (iy == 2) loc_z = curr_loc;

			for( G4int ix = 0; ix < iDims; ix++ ) {
				G4int orient_copy = ix + (iy)*iDims + (iz)*iDims*iDims;
				G4double curr_orient = pp->det_orient_arr[orient_copy];
				if (iy == 0 && ix == 0) u_x = curr_orient;
				else if (iy == 0 && ix == 1) v_x = curr_orient;
				else if (iy == 0 && ix == 2) w_x = curr_orient;
				else if (iy == 1 && ix == 0) u_y = curr_orient;
				else if (iy == 1 && ix == 1) v_y = curr_orient;
				else if (iy == 1 && ix == 2) w_y = curr_orient;
				else if (iy == 2 && ix == 0) u_z = curr_orient;
				else if (iy == 2 && ix == 1) v_z = curr_orient;
				else if (iy == 2 && ix == 2) w_z = curr_orient;
			} // end ix loop
		} // end iy loop
		// creating the center vector of current detector and 3 orientation vectors
		G4ThreeVector loc = G4ThreeVector(loc_x, loc_y, loc_z)*og;
		G4ThreeVector u = G4ThreeVector(u_x, u_y, u_z);
		G4ThreeVector v = G4ThreeVector(v_x, v_y, v_z);
		G4ThreeVector w = G4ThreeVector(w_x, w_y, w_z);
				// create rotation matrix
				G4RotationMatrix rotm  = G4RotationMatrix(u, v, w);
				G4Transform3D transform = G4Transform3D(rotm,loc);
				if (pp->int_map["i_BUILD_DETECTORS"] == 1){
					new G4PVPlacement(transform, //position
									  detectorLV, // logical volume
									  "detector", // name
									  worldLV, //mother volume
									  false, //no boolean operation
									  iz, // copy number
									  false); // check overlaps
		}
	} // end iz loop


	//Pixel
	//G4Box* detectorPixelS = new G4Box("detectorCell",detector_sizeX, detector_sizeY/pp->int_map["i_NUM_OF_DET_COLS"], detector_sizeZ);
	//detectorPixelLV = new G4LogicalVolume(detectorPixelS, detectorMat,"detectorPixel");
	//new G4PVReplica("detectorPixelP",detectorPixelLV,detectorLV,kYAxis,pp->int_map["i_NUM_OF_DET_COLS"],2*detector_sizeY/pp->int_map["i_NUM_OF_DET_COLS"]);


	//world
	//setting visualization attributes to logical elements
	G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	visAttributes->SetVisibility(false);
	worldLV->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes);

	//detector
	visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));
	detectorLV->SetVisAttributes(visAttributes);
	fVisAttributes.push_back(visAttributes);

	//building phantom
	if (pp->int_map["i_BUILD_PHANTOM"] == 1){
		//initialize materials
		if(pp->int_map["i_CT_PHANTOM"]){
			//InitialisationOfMaterialsCT();
			//basic and wrong representation of GT materials
			InitialisationOfMaterialsCT_basic();
			ReadPhantomDataCT();
		}
		else if (pp->int_map["i_XCAT"]){
			//initialize materials
		    InitialisationOfMaterials_XCAT();
		    InitialisationOfMaterialsMap_XCAT();
			ReadPhantomData_XCAT();
		}
		else{
			InitialisationOfMaterials();
			ReadPhantomData();
		}
		ConstructPhantomContainer();
		ConstructPhantom();
	}


	// User Limits

  // Set additional contraints on the track, with G4UserSpecialCuts
  //TODO: Why +5keV?
   G4double maxStep=DBL_MAX, maxLength = DBL_MAX, maxTime = DBL_MAX, minEkin = pp->double_map["f_PARTICLE_ENERGY"]*keV + 5*keV;
   detectorLV->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,minEkin));
//   if (BUILD_PHANTOM == 1){
//	   fvoxel_logic->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,minEkin));
//   }

 //-------------------regions-------------------

//   G4Region* Container_region = new G4Region("ContainerRegion");
//   fContainer_logic->SetRegion(Container_region);
//   Container_region->AddRootLogicalVolume(fContainer_logic);



	//always return the physical World
	return worldPHS;
}

void B1DetectorConstruction::setContainerRotation(G4double delta){
	//sets the rotation of the phantom
	G4RotationMatrix* rot = new G4RotationMatrix();
	rot->rotateZ(delta);
	fContainer_phys->SetRotation(rot);
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
	return;
}

void B1DetectorConstruction::ConstructPhantomContainer()
{
    // read the phantom offset and orientation
	G4double loc_x, loc_y, loc_z, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z;

	for (G4int row=0; row<4; row++){
		for (G4int col=0; col<3; col++){
			G4int curr_idx = row*3 + col;
			G4double curr_val = pp->phantom_loc_orient[curr_idx];
			if (row == 0 && col == 0) loc_x = curr_val;
			else if (row == 0 && col == 1) loc_y = curr_val;
			else if (row == 0 && col == 2) loc_z = curr_val;
			else if (row == 1 && col == 0) u_x = curr_val;
			else if (row == 1 && col == 1) u_y = curr_val;
			else if (row == 1 && col == 2) u_z = curr_val;
			else if (row == 2 && col == 0) v_x = curr_val;
			else if (row == 2 && col == 1) v_y = curr_val;
			else if (row == 2 && col == 2) v_z = curr_val;
			else if (row == 1 && col == 0) w_x = curr_val;
			else if (row == 1 && col == 1) w_y = curr_val;
			else if (row == 1 && col == 2) w_z = curr_val;
		}
	}
	G4double og = 1*mm;
	phantom_loc = G4ThreeVector(loc_x, loc_y, loc_z)*og;
	G4ThreeVector u = G4ThreeVector(u_x, u_y, u_z);
	G4ThreeVector v = G4ThreeVector(v_x, v_y, v_z);
	G4ThreeVector w = G4ThreeVector(w_x, w_y, w_z);

	// create rotation matrix
	phantom_rotm  = G4RotationMatrix(u, v, w);
//	phantom_rotm  = G4RotationMatrix();
//	G4Transform3D transform = G4Transform3D(phantom_rotm,phantom_loc);

  //----- Define the volume that contains all the voxels
  fContainer_solid = new G4Box("phantomContainer",fNVoxelX*fVoxelHalfDimX,
                               fNVoxelY*fVoxelHalfDimY,
                               fNVoxelZ*fVoxelHalfDimZ);
  fContainer_logic =
    new G4LogicalVolume( fContainer_solid,
   //the material is not important, it will be fully filled by the voxels
                         fMaterials[0],
                         "phantomContainer",
                         0, 0, 0 );
  //--- Place it on the world
  G4int copyNum  = 1;
  fContainer_phys =
    new G4PVPlacement(0,
    				  phantom_loc,// orientation//  location
                      fContainer_logic,     // The logic volume
                      "phantomContainer",  // Name
                      worldLV,  // Mother
                      false,           // No op. bool.
                      1);              // Copy number

  //fContainer_logic->SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.,0.)));
}

void B1DetectorConstruction::InitialisationOfMaterialsCT(){
	//initialization based on the materials from DICOM example
	//TODO: remove hard coding!
	G4String MaterialName[10] = {
			"G4_AIR",
			"G4_LUNG_ICRP",
			"G4_ADIPOSE_TISSUE_ICRP",
			"G4_WATER",
			"G4_MUSCLE_WITH_SUCROSE",
			"G4_B-100_BONE",
			"G4_BONE_COMPACT_ICRU",
			"G4_BONE_CORTICAL_ICRP",
			"G4_Fe"
	};
	//const G4MaterialTable* matTab = G4Material::GetMaterialTable();
	//iterate over the base materials

	//TODO: delete later
	std::ofstream output;
	std::string fileName = std::string(pp->string_map["s_OUTPUT_DIR"]) + "/Dicom_base_materials.csv";
	output.open(fileName.c_str());

	for (G4int i = 0; i<pp->int_map["i_NUM_OF_BASE_MATERIALS"]; i++){
		G4Material* material = 0;
		material = G4NistManager::Instance()->FindOrBuildMaterial(MaterialName[i]);
		thePhantomMaterialsOriginal[i] = material;


		const G4ElementVector* curr_element_vector = material->GetElementVector();
		const G4double* curr_frac_vector = material->GetFractionVector();
		G4int nElements = material->GetNumberOfElements();
		output << material->GetDensity()/(g/cm3) << ',';
		for (G4int el=0 ; el<nElements ; el++) {
			G4Element* el_i =  (*curr_element_vector)[el];
			G4String n = el_i->GetName();
			G4double ZZ = el_i->GetZ();
			G4double frac = curr_frac_vector[el];
			output << ZZ << ',' << frac << ',';
		}
		output << '\n';

	}
	output.close();

}

void B1DetectorConstruction::InitialisationOfMaterialsCT_basic()
{
    // Creating elements :
    G4double z, a, density;
    G4int numberofElements;
    G4String name, symbol;

    G4int numOfEl = pp->int_map["i_NUM_OF_ELEMENTS"];

    G4Element** Elements = new G4Element*[numOfEl];
    G4String* ElementName = new G4String[pp->int_map["i_NUM_OF_ELEMENTS"]];

    for (G4int i=0; i<numOfEl; i++){
    	ElementName[i] = "Element" + IntToString(i);
    	a = pp->A_list[i] * g/mole;
		// Elements[i] = G4NistManager::Instance()->FindOrBuildElement(pp->elements[i], false);

    	Elements[i] = new G4Element(name= ElementName[i], symbol = pp->elements[i],
    								pp->Z_list[i], a);
    }

/*
    G4Element* elH = new G4Element( name = "Hydrogen",
                                       symbol = "H",
                                       z = 1.0, a = 1.008  * g/mole );

    G4Element* elC = new G4Element( name = "Carbon",
   								   symbol = "C",
   								   z = 6.0, a = 12.011 * g/mole );

    G4Element* elN = new G4Element( name = "Nitrogen",
                                       symbol = "N",
                                       z = 7.0, a = 14.007 * g/mole );

    G4Element* elO = new G4Element( name = "Oxygen",
                                   symbol = "O",
                                   z = 8.0, a = 16.00  * g/mole );

    G4Element* elP = new G4Element( name = "Phosphorus",
   									symbol = "P",
   									z= 15.0, a = 30.97376 * g/mole );

    G4Element* elCa = new G4Element( name="Calcium",
                                    symbol = "Ca",
                                    z = 20.0, a = 40.078* g/mole );
*/

//*******************************************************************************************************************************************************
    G4String fname;

    G4double* materials_arr = pp->materials_arr;

	//pointers to materials
	G4Material* material;
	G4String* MaterialName = new G4String[pp->int_map["i_NUM_OF_BASE_MATERIALS"]];
	//iterate over the base materials
	G4int mat = 0;

	for (mat = 0; mat<pp->int_map["i_NUM_OF_BASE_MATERIALS"]; mat ++){
		// material - place in an array
		MaterialName[mat] = "mat" + IntToString(mat);
		G4double* fracs = new G4double[numOfEl];
		// read density
		//fin >> density;
		//if (mat==0) {density=1;}
		//else {density = 0.000000000000000001;}
		material = new G4Material( name = MaterialName[mat],
											   1*g/cm3,
											   numberofElements = numOfEl);
//		std::cout << "material number: " << mat << " rho: " << density << " fractions: ";
		// read fractions
		for (G4int i=0; i<numOfEl; i++){
			fracs[i] = materials_arr[numOfEl*mat + i];
			std::cout << fracs[i] << ",";
			material->AddElement(Elements[i],fracs[i]);
		}
		std::cout << "." << std::endl;


		//add material to fMaterials
		thePhantomMaterialsOriginal[mat] = material;
		delete [] fracs;
	}
	delete [] MaterialName;
	delete [] ElementName;
}

void B1DetectorConstruction::InitialisationOfMaterials()
{
    // Creating elements :
    G4double z, a, density;
    G4String name, symbol;
    G4int numberofElements;
    G4Element* elH = new G4Element( name = "Hydrogen",
                                           symbol = "H",
                                           z = 1.0, a = 1.008  * g/mole );

        G4Element* elC = new G4Element( name = "Carbon",
       								   symbol = "C",
       								   z = 6.0, a = 12.011 * g/mole );

        G4Element* elN = new G4Element( name = "Nitrogen",
                                           symbol = "N",
                                           z = 7.0, a = 14.007 * g/mole );

        G4Element* elO = new G4Element( name = "Oxygen",
                                       symbol = "O",
                                       z = 8.0, a = 16.00  * g/mole );

        G4Element* elP = new G4Element( name = "Phosphorus",
       									symbol = "P",
       									z= 15.0, a = 30.97376 * g/mole );

        G4Element* elCa = new G4Element( name="Calcium",
                                        symbol = "Ca",
                                        z = 20.0, a = 40.078* g/mole );

//*******************************************************************************************************************************************************
	G4double* materials_arr = pp->materials_arr;
	//pointers to materials
	G4Material* material;
	G4String* MaterialName = new G4String[pp->int_map["i_NUM_OF_BASE_MATERIALS"]];
	//iterate over the base materials
	G4int mat = 0;
	for (mat = 0; mat<pp->int_map["i_NUM_OF_BASE_MATERIALS"]; mat ++){
		// material - place in an array
		MaterialName[mat] = "mat" + IntToString(mat);
		G4double* fracs = new G4double[pp->int_map["i_NUM_OF_ELEMENTS"]];
		// read density

		material = new G4Material( name = MaterialName[mat],
											   density = density*g/cm3,
											   numberofElements = pp->int_map["i_NUM_OF_ELEMENTS"]);
		std::cout << "material number: " << mat << " rho: " << 1 << " fractions: ";
		// read fractions
		for (G4int i=0; i<pp->int_map["i_NUM_OF_ELEMENTS"]; i++){
			fracs[i] = materials_arr[pp->int_map["i_NUM_OF_ELEMENTS"]*mat + i];
			std::cout << fracs[i] << ",";
		}
		std::cout << "." << std::endl;

		//adding elements according to fractions
				material->AddElement(elH,fracs[0]);
				material->AddElement(elC,fracs[1]);
				material->AddElement(elN,fracs[2]);
				material->AddElement(elO,fracs[3]);
				material->AddElement(elP,fracs[4]);
				material->AddElement(elCa,fracs[5]);

		//add material to fMaterials
		fBaseMaterials[mat] = material;
		delete [] fracs;
	}
	delete [] MaterialName;
	delete [] fMateIDs;
}


void B1DetectorConstruction::ReadPhantomData()
{
    // initiallize fMateIDs
	fMateIDs = new size_t[pp->int_map["i_NUM_OF_VOXELS"]];

    G4String fname = pp->string_map["s_FILE_VOXEL_TO_MATERIALS"];
	std::ifstream fin(fname.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   G4Exception("Can't read voxels_to_materials_file",
					"",
					FatalErrorInArgument,
					G4String("File not found " + fname ).c_str());
	}

    G4double index;
    G4double newRho;
	G4Material* newMaterial;
	G4String baseMaterialName;
	G4String newMaterialName;

	for( G4int iz = 0; iz < pp->int_map["i_NUM_OF_Z_SLICES"]; iz++ ) {
		for( G4int iy = 0; iy < pp->int_map["i_NUM_OF_VOXELS_Y"]; iy++ ) {
			for( G4int ix = 0; ix < pp->int_map["i_NUM_OF_VOXELS_X"]; ix++ ) {
				if( fin.eof() ) break;
				G4int voxel = ix + (iy)*pp->int_map["i_NUM_OF_VOXELS_X"] + (iz)*pp->int_map["i_NUM_OF_VOXELS_X"]*pp->int_map["i_NUM_OF_VOXELS_Y"];
				fin >> index;
				fin >> newRho;
				baseMaterialName = fBaseMaterials[G4int(index)]->GetName();
				//const G4double* fracs = fBaseMaterials[G4int(index)]->GetFractionVector();
				//std::cout << "voxel number: " << voxel << " rho: " << newRho << " atom1 fraction: " << fracs[0] << " atom 2 fraction: " << fracs[1] << "\n";
				newMaterialName = "voxelMat" + IntToString(voxel);
				newMaterial = G4NistManager::Instance()->
				  BuildMaterialWithNewDensity(newMaterialName,baseMaterialName,newRho*g/cm3);

				fMaterials.push_back(newMaterial);
				fMateIDs[voxel] = voxel;
			}
		}
	}
}


void B1DetectorConstruction::ReadPhantomDataCT()
{

  //---- Extract number of voxels and voxel dimensions
  fNVoxelX = pp->int_map["i_NUM_OF_VOXELS_X"];
  fNVoxelY = pp->int_map["i_NUM_OF_VOXELS_Y"];
  fNVoxelZ = pp->int_map["i_NUM_OF_Z_SLICES"];

  fVoxelHalfDimX = pp->double_map["f_VOXEL_HALF_X"]*mm;
  fVoxelHalfDimY = pp->double_map["f_VOXEL_HALF_Y"]*mm;
  fVoxelHalfDimZ = pp->double_map["f_VOXEL_HALF_Z"]*mm;

  fMateIDs = new size_t[fNVoxelX*fNVoxelY*fNVoxelZ];
  for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
    for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
      for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
        G4int mateID;
        G4int nnew = ix + (iy)*fNVoxelX + (iz)*fNVoxelX*fNVoxelY;
        mateID = pp->ids_arr[nnew];
        if( mateID < 0 || mateID >= pp->int_map["i_NUM_OF_BASE_MATERIALS"] ) {
          G4Exception("GmReadPhantomG4Geometry::ReadPhantomData",
                      "Wrong index in phantom file",
                      FatalException,
                      G4String("It should be between 0 and "
                               + G4UIcommand::ConvertToString(pp->int_map["i_NUM_OF_BASE_MATERIALS"]-1)
                               + ", while it is "
                               + G4UIcommand::ConvertToString(mateID)).c_str());
        }
        fMateIDs[nnew] = mateID;
      }
    }
  }
  ReadVoxelDensities( );

}


void B1DetectorConstruction::ReadVoxelDensities( )
{

	G4String stemp;
	  std::map<G4int, std::pair<G4double,G4double> > densiMinMax;
	  std::map<G4int, std::pair<G4double,G4double> >::iterator mpite;
	  for( size_t ii = 0; ii < thePhantomMaterialsOriginal.size(); ii++ ){
	    densiMinMax[ii] = std::pair<G4double,G4double>(DBL_MAX,-DBL_MAX);
	  }

	  //char* part = getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
	  //G4double densityDiff = -1.;
	  G4double densityDiff = 0.01;
	  //if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);

	  std::map<G4int,G4double> densityDiffs;
	  for( size_t ii = 0; ii < thePhantomMaterialsOriginal.size(); ii++ ){
	    densityDiffs[ii] = densityDiff; //currently all materials with same step
	  }
	  //  densityDiffs[0] = 0.0001; //air


	  //--- Calculate the average material density for each material/density bin
	  std::map< std::pair<G4Material*,G4int>, matInfo* > newMateDens;

	  //---- Read the material densities
	  G4double dens;
	  for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
	    for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
	      for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
	        G4int copyNo = ix + (iy)*fNVoxelX + (iz)*fNVoxelX*fNVoxelY;
	    	dens = pp->dens_arr[copyNo];

	        //--- store the minimum and maximum density for each material (just for printing)
	        mpite = densiMinMax.find( fMateIDs[copyNo] );
	        if( dens < (*mpite).second.first ) (*mpite).second.first = dens;
	        if( dens > (*mpite).second.second ) (*mpite).second.second = dens;
	        //--- Get material from original list of material in file
	        int mateID = fMateIDs[copyNo];
	        std::map<G4int,G4Material*>::const_iterator imite =
	         thePhantomMaterialsOriginal.find(mateID);
	        //        G4cout << copyNo << " mateID " << mateID << G4endl;
	        //--- Check if density is equal to the original material density
			G4double curr_dens = (*imite).second->GetDensity()/CLHEP::g*CLHEP::cm3; 
	        if( std::fabs(dens - curr_dens) < 1.e-9 ) continue;

	        //--- Build material name with thePhantomMaterialsOriginal name + density
	        //        float densityBin = densityDiffs[mateID] * (G4int(dens/densityDiffs[mateID])+0.5);
	        G4int densityBin = (G4int(dens/densityDiffs[mateID]));

	        G4String mateName = (*imite).second->GetName()+G4UIcommand::ConvertToString(densityBin);
	        //--- Look if it is the first voxel with this material/densityBin
	        std::pair<G4Material*,G4int> matdens((*imite).second, densityBin );

	        std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite =
	         newMateDens.find( matdens );
	        if( mppite != newMateDens.end() ){
	          matInfo* mi = (*mppite).second;
	          mi->fSumdens += dens;
	          mi->fNvoxels++;
	          fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
	        } else {
	          matInfo* mi = new matInfo;
	          mi->fSumdens = dens;
	          mi->fNvoxels = 1;
	          mi->fId = newMateDens.size()+1;
	          newMateDens[matdens] = mi;
	          fMateIDs[copyNo] = thePhantomMaterialsOriginal.size()-1 + mi->fId;
	        }
	      }
	    }
	  }

	  if( densityDiff != -1. ) {
	    for( mpite = densiMinMax.begin(); mpite != densiMinMax.end(); mpite++ ){

	    }
	  }

	  //----- Build the list of phantom materials that go to Parameterisation
	  //--- Add original materials
	  std::map<G4int,G4Material*>::const_iterator mimite;
	  for( mimite = thePhantomMaterialsOriginal.begin(); mimite != thePhantomMaterialsOriginal.end();
	   mimite++ ){
	    fMaterials.push_back( (*mimite).second );
	  }
	  //
	  //---- Build and add new materials

	  //G4int sDens = newMateDens.sisze();
	  std::map< G4int, std::pair<G4Material*,G4int> > sorted_map_mat;
	  std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite;

	  for( mppite= newMateDens.begin(); mppite != newMateDens.end(); mppite++ ){
		  G4int idd = (*mppite).second->fId;
		  sorted_map_mat[idd] = (*mppite).first;
	  }

	  std::map< G4int, std::pair<G4Material*,G4int> >::iterator iterator_sorted;
	  for( iterator_sorted= sorted_map_mat.begin(); iterator_sorted != sorted_map_mat.end(); iterator_sorted++ ){
		  std::pair<G4Material*,G4int> key = (*iterator_sorted).second;
		  mppite = newMateDens.find( key );

		  G4double averdens = (*mppite).second->fSumdens/(*mppite).second->fNvoxels;
		  G4double saverdens = G4int(1000.001*averdens)/1000.;

		  G4String mateName = ((*mppite).first).first->GetName() + "_"
		      + G4UIcommand::ConvertToString(saverdens);
		  fMaterials.push_back( BuildMaterialWithChangingDensity(
		      (*mppite).first.first, averdens, mateName ) );
	  }

/*
	  for( mppite= newMateDens.begin(); mppite != newMateDens.end(); mppite++ ){
	    G4double averdens = (*mppite).second->fSumdens/(*mppite).second->fNvoxels;
	    G4double saverdens = G4int(1000.001*averdens)/1000.;

	      G4String mateName = ((*mppite).first).first->GetName() + "_"
	       + G4UIcommand::ConvertToString(saverdens);
	    fMaterials.push_back( BuildMaterialWithChangingDensity(
	     (*mppite).first.first, averdens, mateName ) );
	  }
*/
}

G4Material* B1DetectorConstruction::BuildMaterialWithChangingDensity(
           const G4Material* origMate, float density, G4String newMateName )
{
  //----- Copy original material, but with new density
  G4int nelem = origMate->GetNumberOfElements();
  G4Material* mate = new G4Material( newMateName, density*g/cm3, nelem,
                                     kStateUndefined, STP_Temperature );

  for( G4int ii = 0; ii < nelem; ii++ ){
    G4double frac = origMate->GetFractionVector()[ii];
    G4Element* elem = const_cast<G4Element*>(origMate->GetElement(ii));
    mate->AddElement( elem, frac );
  }

  return mate;
}

//**************************************************************************************************************************XCAT
void B1DetectorConstruction::InitialisationOfMaterialsMap_XCAT()
{
	G4String fname = pp->string_map["s_FILE_XCAT_ID_TO_COMP"];
	std::ifstream fin(fname.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   G4Exception("B1DetectorConstruction::ReadPhantomDataFile",
	                "",
	                FatalErrorInArgument,
	                G4String("File not found " + fname ).c_str());
	  }
	//read header - dummy line
	std::string line;
	std::getline(fin,line);
	//std::cout << line << std::endl;
	//read file and create a map
	while (true) {
	    if( fin.eof() ) break;

	    G4int NumOfMaterials;
	    G4int ID;
	    G4int MateA;
	    G4int MateB;

	    fin >> NumOfMaterials;
	    fin >> ID;
	    //std::cout<<"NumOfMaterials: " << NumOfMaterials << " ID: " << ID << std::endl;
	    //create material struct
	    materialIDs MaterialStruct;
	    if (NumOfMaterials==1){
	    	fin >> MateA;
	    	//std::cout<<"MateA: " << MateA << std::endl;
	    	MaterialStruct.mat1ID=MateA;
	    	MaterialStruct.mat2ID=-1;
	    } else { //NumOfMaterials==2
	    	fin >> MateA;
	    	fin >> MateB;
	    	//std::cout<<"MateA: " << MateA << "MateB: " << MateB << std::endl;
	    	MaterialStruct.mat1ID=MateA;
	    	MaterialStruct.mat2ID=MateB;
	    }

	    //insert to map
	    std::pair<std::map<G4int,materialIDs>::iterator,bool> ret;
	    //std::cout << "inserting ID: " << ID << std::endl;
	    ret = intensityToMateID.insert ( std::pair<G4int,materialIDs>(ID,MaterialStruct) );
	    if (ret.second==false) {
	      std::cout << "element " << ID << " already existed" << std::endl;
	      //std::cout << " with a value of " << intensityToMateID[ret.first] << '\n';
	    }
	}
	//std::cout << "intensityToMateID.size() is " << intensityToMateID.size() << std::endl;

}

void B1DetectorConstruction::InitialisationOfMaterials_XCAT()
{
    // Creating elements :
    G4double z, a, density;
    G4String name, symbol;


    G4Element* elH = new G4Element( name = "Hydrogen",
                                   symbol = "H",
                                   z = 1.0, a = 1.008  * g/mole );
    G4Element* elHe = new G4Element( name = "Helium",
                                   symbol = "He",
                                   z = 2.0, a = 4.0026  * g/mole );
    G4Element* elLi = new G4Element( name = "Lithium",
                                   symbol = "Li",
                                   z = 3.0, a = 6.941  * g/mole );
    G4Element* elBe = new G4Element( name = "Beryllium",
                                   symbol = "Be",
                                   z = 4.0, a = 9.012182  * g/mole );
    G4Element* elB = new G4Element( name = "Boron",
								   symbol = "B",
								   z = 5.0, a = 10.811  * g/mole );
    G4Element* elC = new G4Element( name = "Carbon",
								   symbol = "C",
								   z = 6.0, a = 12.011 * g/mole );
    G4Element* elN = new G4Element( name = "Nitrogen",
                                   symbol = "N",
                                   z = 7.0, a = 14.007 * g/mole );
    G4Element* elO = new G4Element( name = "Oxygen",
                                   symbol = "O",
                                   z = 8.0, a = 16.00  * g/mole );
    G4Element* elF = new G4Element( name = "Fluorine",
								   symbol = "F",
								   z = 9.0, a = 18.998404  * g/mole );
    G4Element* elNe = new G4Element( name = "Neon",
								   symbol = "Ne",
								   z = 10.0, a = 20.1797  * g/mole );
    G4Element* elNa = new G4Element( name = "Sodium",
                                    symbol = "Na",
                                    z= 11.0, a = 22.98977 * g/mole );
    G4Element* elMg = new G4Element( name = "Magnesium",
									symbol = "Mg",
									z= 12.0, a = 24.305 * g/mole );
    G4Element* elAl = new G4Element( name = "Aluminum",
									symbol = "Al",
									z= 13.0, a = 26.981539 * g/mole );
    G4Element* elP = new G4Element( name = "Phosphorus",
									symbol = "P",
									z= 15.0, a = 30.97376 * g/mole );
    G4Element* elS = new G4Element( name = "Sulfur",
                                   	symbol = "S",
                                   	z = 16.0,a = 32.065* g/mole );
    G4Element* elCl = new G4Element( name = "Chlorine",
                                    symbol = "Cl",
                                    z = 17.0, a = 35.453* g/mole );
    G4Element* elAr = new G4Element( name = "Argon",
									symbol = "Ar",
									z= 18.0, a = 39.948 * g/mole );
    G4Element* elK = new G4Element( name = "Potassium",
                                   	symbol = "K",
                                   	z = 19.0, a = 39.0983* g/mole );
    G4Element* elCa = new G4Element( name="Calcium",
                                    symbol = "Ca",
                                    z = 20.0, a = 40.078* g/mole );
    G4Element* elSc = new G4Element( name="Scandium",
                                    symbol = "Sc",
                                    z = 21.0, a = 44.95591 * g/mole );
    G4Element* elTi = new G4Element( name="Titanium",
                                    symbol = "Ti",
                                    z = 22.0, a = 47.867 * g/mole );
    G4Element* elV = new G4Element( name="Vanadium",
                                    symbol = "V",
                                    z = 23.0, a = 50.9415 * g/mole );
    G4Element* elCr = new G4Element( name="Chromium",
                                    symbol = "Cr",
                                    z = 24.0, a = 51.9961 * g/mole );
    G4Element* elMn = new G4Element( name="Manganese",
                                    symbol = "Mn",
                                    z = 25.0, a = 54.93805 * g/mole );
    G4Element* elFe = new G4Element( name = "Iron",
                                    symbol = "Fe",
                                    z = 26, a = 55.845* g/mole );
    G4Element* elI = new G4Element( name = "Iodine",
                                    symbol = "I",
                                    z = 53, a = 126.90447 * g/mole );
    G4Element* elPb = new G4Element( name = "Lead",
                                    symbol = "Pb",
                                    z = 82, a = 207.2 * g/mole );


    // Creating Materials :
    G4int numberofElements;

    // Water
    G4Material* water = new G4Material( "Water",
                                       density = 1.0*g/cm3,
                                       numberofElements = 2 );
    water->AddElement(elH,0.112);
    water->AddElement(elO,0.888);

    // Muscle
	G4Material* muscle = new G4Material( "Muscle",
										density = 1.05*g/cm3,
										numberofElements = 9 );
	muscle->AddElement(elH,0.102);
	muscle->AddElement(elC,0.143);
	muscle->AddElement(elN,0.034);
	muscle->AddElement(elO,0.710);
	muscle->AddElement(elNa,0.001);
	muscle->AddElement(elP,0.002);
	muscle->AddElement(elS,0.003);
	muscle->AddElement(elCl,0.001);
	muscle->AddElement(elK,0.004);

    //  Lung Inhale
    G4Material* lung = new G4Material( "Lung",
                                            density = 0.29*g/cm3,
                                            numberofElements = 9);
    lung->AddElement(elH,0.103);
    lung->AddElement(elC,0.105);
    lung->AddElement(elN,0.031);
    lung->AddElement(elO,0.749);
    lung->AddElement(elNa,0.002);
    lung->AddElement(elP,0.002);
    lung->AddElement(elS,0.003);
    lung->AddElement(elCl,0.003);
    lung->AddElement(elK,0.002);

    // Dry spine
    G4Material* dry_spine = new G4Material( "DrySpine",
                                            density = 1.42*g/cm3,
                                            numberofElements = 11);
    dry_spine->AddElement(elH,0.063);
    dry_spine->AddElement(elC,0.261);
    dry_spine->AddElement(elN,0.039);
    dry_spine->AddElement(elO,0.436);
    dry_spine->AddElement(elNa,0.001);
    dry_spine->AddElement(elMg,0.001);
    dry_spine->AddElement(elP,0.061);
    dry_spine->AddElement(elS,0.003);
    dry_spine->AddElement(elCl,0.001);
    dry_spine->AddElement(elK,0.001);
    dry_spine->AddElement(elCa,0.133);

    // Dry rib
    G4Material* dry_rib = new G4Material( "DryRib",
                                            density = 1.92*g/cm3,
                                            numberofElements = 9);
    dry_rib->AddElement(elH,0.034);
    dry_rib->AddElement(elC,0.155);
    dry_rib->AddElement(elN,0.042);
    dry_rib->AddElement(elO,0.435);
    dry_rib->AddElement(elNa,0.001);
    dry_rib->AddElement(elMg,0.002);
    dry_rib->AddElement(elP,0.103);
    dry_rib->AddElement(elS,0.003);
    dry_rib->AddElement(elCa,0.225);


    // Adipose tissue
    G4Material* adiposeTissue = new G4Material( "AdiposeTissue",
                                               density = 0.95*g/cm3,
                                               numberofElements = 7);
    adiposeTissue->AddElement(elH,0.114);
    adiposeTissue->AddElement(elC,0.598);
    adiposeTissue->AddElement(elN,0.007);
    adiposeTissue->AddElement(elO,0.278);
    adiposeTissue->AddElement(elNa,0.001);
    adiposeTissue->AddElement(elS,0.001);
    adiposeTissue->AddElement(elCl,0.001);

    // Blood
    G4Material* blood = new G4Material( "Blood",
                                               density = 1.06*g/cm3,
                                               numberofElements = 10);
    blood->AddElement(elH,0.102);
    blood->AddElement(elC,0.11);
    blood->AddElement(elN,0.033);
    blood->AddElement(elO,0.745);
    blood->AddElement(elNa,0.001);
    blood->AddElement(elP,0.001);
    blood->AddElement(elS,0.002);
    blood->AddElement(elCl,0.003);
    blood->AddElement(elK,0.002);
    blood->AddElement(elFe,0.001);

    // Heart
    G4Material* heart = new G4Material( "Heart",
                                               density = 1.05*g/cm3,
                                               numberofElements = 9);
    heart->AddElement(elH,0.104);
    heart->AddElement(elC,0.139);
    heart->AddElement(elN,0.029);
    heart->AddElement(elO,0.718);
    heart->AddElement(elNa,0.001);
    heart->AddElement(elP,0.002);
    heart->AddElement(elS,0.002);
    heart->AddElement(elCl,0.002);
    heart->AddElement(elK,0.003);

    // Kidney
    G4Material* kidney = new G4Material( "Kidney",
                                               density = 1.05*g/cm3,
                                               numberofElements = 10);
    kidney->AddElement(elH,0.103);
    kidney->AddElement(elC,0.132);
    kidney->AddElement(elN,0.03);
    kidney->AddElement(elO,0.724);
    kidney->AddElement(elNa,0.002);
    kidney->AddElement(elP,0.002);
    kidney->AddElement(elS,0.002);
    kidney->AddElement(elCl,0.002);
    kidney->AddElement(elK,0.002);
    kidney->AddElement(elCa,0.001);

    // Liver
	G4Material* liver = new G4Material( "Liver",
									   density = 1.06*g/cm3,
									   numberofElements = 9);
	liver->AddElement(elH,0.102);
	liver->AddElement(elC,0.139);
	liver->AddElement(elN,0.030);
	liver->AddElement(elO,0.716);
	liver->AddElement(elNa,0.002);
	liver->AddElement(elP,0.003);
	liver->AddElement(elS,0.003);
	liver->AddElement(elCl,0.002);
	liver->AddElement(elK,0.003);

    // Lymph
	G4Material* lymph = new G4Material( "Lymph",
									   density = 1.03*g/cm3,
									   numberofElements = 7);
	lymph->AddElement(elH,0.108);
	lymph->AddElement(elC,0.041);
	lymph->AddElement(elN,0.011);
	lymph->AddElement(elO,0.832);
	lymph->AddElement(elNa,0.003);
	lymph->AddElement(elS,0.001);
	lymph->AddElement(elCl,0.004);

    // Pancreas
	G4Material* pancreas = new G4Material( "Pancreas",
									   density = 1.04*g/cm3,
									   numberofElements = 9);
	pancreas->AddElement(elH,0.106);
	pancreas->AddElement(elC,0.169);
	pancreas->AddElement(elN,0.022);
	pancreas->AddElement(elO,0.694);
	pancreas->AddElement(elNa,0.002);
	pancreas->AddElement(elP,0.002);
	pancreas->AddElement(elS,0.001);
	pancreas->AddElement(elCl,0.002);
	pancreas->AddElement(elK,0.002);

    // Intestine
	G4Material* intestine = new G4Material( "Intestine",
									   density = 1.03*g/cm3,
									   numberofElements = 9);
	intestine->AddElement(elH,0.106);
	intestine->AddElement(elC,0.115);
	intestine->AddElement(elN,0.022);
	intestine->AddElement(elO,0.751);
	intestine->AddElement(elNa,0.001);
	intestine->AddElement(elP,0.001);
	intestine->AddElement(elS,0.001);
	intestine->AddElement(elCl,0.002);
	intestine->AddElement(elK,0.001);

    // Skull
	G4Material* skull = new G4Material( "Skull",
									   density = 1.61*g/cm3,
									   numberofElements = 9);
	skull->AddElement(elH,0.05);
	skull->AddElement(elC,0.212);
	skull->AddElement(elN,0.04);
	skull->AddElement(elO,0.435);
	skull->AddElement(elNa,0.001);
	skull->AddElement(elMg,0.002);
	skull->AddElement(elP,0.081);
	skull->AddElement(elS,0.003);
	skull->AddElement(elCa,0.176);

    // Cartilage
	G4Material* cartilage = new G4Material( "Cartilage",
									   density = 1.10*g/cm3,
									   numberofElements = 8);
	cartilage->AddElement(elH,0.096);
	cartilage->AddElement(elC,0.099);
	cartilage->AddElement(elN,0.022);
	cartilage->AddElement(elO,0.744);
	cartilage->AddElement(elNa,0.005);
	cartilage->AddElement(elP,0.022);
	cartilage->AddElement(elS,0.009);
	cartilage->AddElement(elCl,0.003);

    // Brain
	G4Material* brain = new G4Material( "Brain",
									   density = 1.04*g/cm3,
									   numberofElements = 9);
	brain->AddElement(elH,0.107);
	brain->AddElement(elC,0.145);
	brain->AddElement(elN,0.022);
	brain->AddElement(elO,0.712);
	brain->AddElement(elNa,0.002);
	brain->AddElement(elP,0.004);
	brain->AddElement(elS,0.002);
	brain->AddElement(elCl,0.003);
	brain->AddElement(elK,0.003);

    // Spleen
	G4Material* spleen = new G4Material( "Spleen",
									   density = 1.06*g/cm3,
									   numberofElements = 9);
	spleen->AddElement(elH,0.103);
	spleen->AddElement(elC,0.113);
	spleen->AddElement(elN,0.032);
	spleen->AddElement(elO,0.741);
	spleen->AddElement(elNa,0.001);
	spleen->AddElement(elP,0.003);
	spleen->AddElement(elS,0.002);
	spleen->AddElement(elCl,0.002);
	spleen->AddElement(elK,0.003);

    // Iodine Blood
	G4Material* iodine_blood = new G4Material( "IodineBlood",
									   density = 1.09096*g/cm3,
									   numberofElements = 11);
	iodine_blood->AddElement(elH,0.101184);
	iodine_blood->AddElement(elC,0.10912);
	iodine_blood->AddElement(elN,0.032736);
	iodine_blood->AddElement(elO,0.73904);
	iodine_blood->AddElement(elNa,0.000992);
	iodine_blood->AddElement(elP,0.000992);
	iodine_blood->AddElement(elS,0.001984);
	iodine_blood->AddElement(elCl,0.002976);
	iodine_blood->AddElement(elK,0.001984);
	iodine_blood->AddElement(elFe,0.000992);
	iodine_blood->AddElement(elI,0.008);

    // Iron
    G4Material* iron = new G4Material( "Iron",
                                        density = 7.86 * g/cm3,
                                        numberofElements = 1 );
    iron->AddElement(elFe,1.0);

    // Pmma
    G4Material* pmma = new G4Material( "Pmma",
                                        density = 1.19 * g/cm3,
                                        numberofElements = 3 );
    pmma->AddElement(elH,0.08);
    pmma->AddElement(elC,0.6);
    pmma->AddElement(elO,0.32);

    // Aluminum
    G4Material* aluminum = new G4Material( "Aluminum",
                                        density = 2.6941 * g/cm3,
                                        numberofElements = 1 );
    aluminum->AddElement(elAl,1.0);

    // Titanium
    G4Material* titanium = new G4Material( "Titanium",
                                        density = 4.53 * g/cm3,
                                        numberofElements = 1 );
    titanium->AddElement(elTi,1.0);

    // Air
    // Get nist Air instead of phantoms air
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
//    G4Material* air = new G4Material( "Air",
//    								0.001205*mg/cm3,
//    								numberofElements = 4 );
//    air->AddElement(elC, 0.000124);
//    air->AddElement(elN, 0.755268);
//    air->AddElement(elO, 0.231781);
//    air->AddElement(elCl, 0.012827);

    // Graphite
    G4Material* graphite = new G4Material( "Graphite",
                                        density = 1.82 * g/cm3,
                                        numberofElements = 1 );
    graphite->AddElement(elC,1.0);

    // Lead
    G4Material* lead = new G4Material( "Lead",
                                        density = 11.34 * g/cm3,
                                        numberofElements = 1 );
    lead->AddElement(elPb,1.0);

    // Breast Mammary
	G4Material* breast_mammary = new G4Material( "breastMammary",
									   density = 1.02*g/cm3,
									   numberofElements = 8);
	breast_mammary->AddElement(elH,0.106);
	breast_mammary->AddElement(elC,0.332);
	breast_mammary->AddElement(elN,0.03);
	breast_mammary->AddElement(elO,0.527);
	breast_mammary->AddElement(elNa,0.001);
	breast_mammary->AddElement(elP,0.001);
	breast_mammary->AddElement(elS,0.002);
	breast_mammary->AddElement(elCl,0.001);

    // Skin
	G4Material* skin = new G4Material( "Skin",
									   density = 1.09*g/cm3,
									   numberofElements = 9);
	skin->AddElement(elH,0.10);
	skin->AddElement(elC,0.204);
	skin->AddElement(elN,0.042);
	skin->AddElement(elO,0.645);
	skin->AddElement(elNa,0.002);
	skin->AddElement(elP,0.001);
	skin->AddElement(elS,0.002);
	skin->AddElement(elCl,0.003);
	skin->AddElement(elK,0.001);

	 // Iodine
	G4Material* iodine = new G4Material( "Iodine",
									   density = 4.933*g/cm3,
									   numberofElements = 1);
	iodine->AddElement(elI,1.0);

	 // eye_lens
	G4Material* eye_lens = new G4Material( "EyeLens",
									   density = 1.07*g/cm3,
									   numberofElements = 8);
	eye_lens->AddElement(elH,0.096);
	eye_lens->AddElement(elC,0.195);
	eye_lens->AddElement(elN,0.057);
	eye_lens->AddElement(elO,0.646);
	eye_lens->AddElement(elNa,0.001);
	eye_lens->AddElement(elP,0.001);
	eye_lens->AddElement(elS,0.003);
	eye_lens->AddElement(elCl,0.001);

	 // ovary
	G4Material* ovary = new G4Material( "Ovary",
									   density = 1.05*g/cm3,
									   numberofElements = 9);
	ovary->AddElement(elH,0.105);
	ovary->AddElement(elC,0.093);
	ovary->AddElement(elN,0.024);
	ovary->AddElement(elO,0.768);
	ovary->AddElement(elNa,0.002);
	ovary->AddElement(elP,0.002);
	ovary->AddElement(elS,0.002);
	ovary->AddElement(elCl,0.002);
	ovary->AddElement(elK,0.002);

	 // red_marrow
	G4Material* red_marrow = new G4Material( "RedMarrow",
									   density = 1.03*g/cm3,
									   numberofElements = 9);
	red_marrow->AddElement(elH,0.105);
	red_marrow->AddElement(elC,0.414);
	red_marrow->AddElement(elN,0.034);
	red_marrow->AddElement(elO,0.439);
	red_marrow->AddElement(elP,0.001);
	red_marrow->AddElement(elS,0.002);
	red_marrow->AddElement(elCl,0.002);
	red_marrow->AddElement(elK,0.002);
	red_marrow->AddElement(elFe,0.001);

	 // yellow_marrow
	G4Material* yellow_marrow = new G4Material( "YellowMarrow",
									   density = 0.98*g/cm3,
									   numberofElements = 7);
	yellow_marrow->AddElement(elH,0.115);
	yellow_marrow->AddElement(elC,0.644);
	yellow_marrow->AddElement(elN,0.007);
	yellow_marrow->AddElement(elO,0.231);
	yellow_marrow->AddElement(elNa,0.001);
	yellow_marrow->AddElement(elS,0.001);
	yellow_marrow->AddElement(elCl,0.001);

	 // testis
	G4Material* testis = new G4Material( "Testis",
									   density = 1.04*g/cm3,
									   numberofElements = 9);
	testis->AddElement(elH,0.106);
	testis->AddElement(elC,0.099);
	testis->AddElement(elN,0.02);
	testis->AddElement(elO,0.766);
	testis->AddElement(elNa,0.002);
	testis->AddElement(elP,0.001);
	testis->AddElement(elS,0.002);
	testis->AddElement(elCl,0.002);
	testis->AddElement(elK,0.002);

	 // thyroid
	G4Material* thyroid = new G4Material( "Thyroid",
									   density = 1.05*g/cm3,
									   numberofElements = 10);
	thyroid->AddElement(elH,0.104);
	thyroid->AddElement(elC,0.119);
	thyroid->AddElement(elN,0.024);
	thyroid->AddElement(elO,0.745);
	thyroid->AddElement(elNa,0.002);
	thyroid->AddElement(elP,0.001);
	thyroid->AddElement(elS,0.001);
	thyroid->AddElement(elCl,0.002);
	thyroid->AddElement(elK,0.001);
	thyroid->AddElement(elI,0.001);

	 // trabecular
	G4Material* trabecular = new G4Material( "trabecular",
									   density = 1.14*g/cm3,
									   numberofElements = 5);
	trabecular->AddElement(elH,0.079);
	trabecular->AddElement(elC,0.6379);
	trabecular->AddElement(elN,0.0423);
	trabecular->AddElement(elO,0.0988);
	trabecular->AddElement(elCa,0.142);

	 // bladder
	G4Material* bladder = new G4Material( "Bladder",
									   density = 1.04*g/cm3,
									   numberofElements = 9);
	bladder->AddElement(elH,0.105);
	bladder->AddElement(elC,0.096);
	bladder->AddElement(elN,0.026);
	bladder->AddElement(elO,0.761);
	bladder->AddElement(elNa,0.002);
	bladder->AddElement(elP,0.002);
	bladder->AddElement(elS,0.002);
	bladder->AddElement(elCl,0.003);
	bladder->AddElement(elK,0.003);

	//TODO: am I calculating the density correctly?
	 // dry_spine with water 0,3
   G4Material* dry_spine_water = new G4Material("dry_spine_water",
                                           density = 1.4*g/cm3,
                                           numberofElements = 9);
   dry_spine_water->AddElement(elH,0.034);
   dry_spine_water->AddElement(elC,0.155);
   dry_spine_water->AddElement(elN,0.042);
   dry_spine_water->AddElement(elO,0.435);
   dry_spine_water->AddElement(elNa,0.001);
   dry_spine_water->AddElement(elMg,0.002);
   dry_spine_water->AddElement(elP,0.103);
   dry_spine_water->AddElement(elS,0.003);
   dry_spine_water->AddElement(elCa,0.225);

	 // dry_rib with water 0,4
	G4Material* dry_rib_water = new G4Material( "dry_rib_water",
												density = 1.55*g/cm3,
												numberofElements = 9);
	dry_rib_water->AddElement(elH,0.034);
	dry_rib_water->AddElement(elC,0.155);
	dry_rib_water->AddElement(elN,0.042);
	dry_rib_water->AddElement(elO,0.435);
	dry_rib_water->AddElement(elNa,0.001);
	dry_rib_water->AddElement(elMg,0.002);
	dry_rib_water->AddElement(elP,0.103);
	dry_rib_water->AddElement(elS,0.003);
	dry_rib_water->AddElement(elCa,0.225);
	 // skull with water 0,13
	G4Material* skull_water = new G4Material( "skull_water",
									   density = (0.5*water->GetDensity())*g/cm3 + (0.5*skull->GetDensity())*g/cm3,
									   numberofElements = 2);
	skull_water->AddMaterial(skull,50.*perCent);
	skull_water->AddMaterial(water,50.*perCent);



    //----- Put the materials in a vector
    fMaterials.push_back(water);             //0
    fMaterials.push_back(muscle);            //1
    fMaterials.push_back(lung);              //2
    fMaterials.push_back(dry_spine);         //3
    fMaterials.push_back(dry_rib);			 //4
    fMaterials.push_back(adiposeTissue);	 //5
    fMaterials.push_back(blood);			 //6
    fMaterials.push_back(heart);			 //7
    fMaterials.push_back(kidney);			 //8
    fMaterials.push_back(liver);			 //9
    fMaterials.push_back(lymph);			 //10
    fMaterials.push_back(pancreas);			 //11
    fMaterials.push_back(intestine);		 //12
    fMaterials.push_back(skull);			 //13
    fMaterials.push_back(cartilage);		 //14
    fMaterials.push_back(brain);			 //15
    fMaterials.push_back(spleen);			 //16
    fMaterials.push_back(iodine_blood);		 //17
    fMaterials.push_back(iron);				 //18
    fMaterials.push_back(pmma);				 //19
    fMaterials.push_back(aluminum);			 //20
    fMaterials.push_back(titanium);			 //21
    fMaterials.push_back(air);				 //22
    fMaterials.push_back(graphite);			 //23
    fMaterials.push_back(lead);				 //24
    fMaterials.push_back(breast_mammary);	 //25
    fMaterials.push_back(skin);				 //26
    fMaterials.push_back(iodine);			 //27
    fMaterials.push_back(eye_lens);			 //28
    fMaterials.push_back(ovary);			 //29
    fMaterials.push_back(red_marrow);		 //30
    fMaterials.push_back(yellow_marrow);	 //31
    fMaterials.push_back(testis);			 //32
    fMaterials.push_back(thyroid);			 //33
    fMaterials.push_back(trabecular);		 //34
    fMaterials.push_back(bladder);			 //35
    fMaterials.push_back(dry_spine_water);	 //36
    fMaterials.push_back(dry_rib_water);	 //37
    fMaterials.push_back(skull_water);		 //38


}

void B1DetectorConstruction::ReadPhantomData_XCAT()
{

	  //---- Extract number of voxels and voxel dimensions
	  fNVoxelX = pp->int_map["i_NUM_OF_VOXELS_X"];
	  fNVoxelY = pp->int_map["i_NUM_OF_VOXELS_Y"];
	  fNVoxelZ = pp->int_map["i_NUM_OF_Z_SLICES"];

	  fVoxelHalfDimX = pp->double_map["f_VOXEL_HALF_X"]*mm;
	  fVoxelHalfDimY = pp->double_map["f_VOXEL_HALF_Y"]*mm;
	  fVoxelHalfDimZ = pp->double_map["f_VOXEL_HALF_Z"]*mm;


    G4int x=10000; //starting value for file names
    //G4int z_shift = 330; //hand
    G4int z_shift = 215; //knee
    for(G4int i = z_shift; i < pp->int_map["i_NUM_OF_Z_SLICES"] + z_shift; i++ ) {
        //--- Read one data file
        G4String fileName = pp->string_map["s_FILE_XCAT_SLICE_PREFIX"] + IntToString(x+i) + ".txt";
        ReadPhantomDataFile_XCAT(fileName,i - z_shift);
    }
}

void B1DetectorConstruction::ReadPhantomDataFile_XCAT(const G4String& fname, G4int sliceNumber)
{

  std::cout << " B1DetectorConstruction::ReadPhantomDataFile opening file "
		 << fname << std::endl;
  //TODO: handle reading from phantom files
  std::ifstream fin(fname.c_str(), std::ios_base::in);
  if( !fin.is_open() ) {
    G4Exception("B1DetectorConstruction::ReadPhantomDataFile",
                "",
                FatalErrorInArgument,
                G4String("File not found " + fname ).c_str());
  }

  //--- If first slice, initiliaze fMateIDs
  if( sliceNumber==0 ) {
    fMateIDs = new size_t[pp->int_map["i_NUM_OF_VOXELS"]];
  }

//comment when no need to get GT
/*
  std::ofstream output;
  std::string fileName = "materails_XCAT_ID" + IntToString(sliceNumber) + ".csv";
  output.open(fileName.c_str());
*/
  G4double mateID;
  // number of voxels from previously read slices
  G4int voxelCopyNo = (sliceNumber)*pp->int_map["i_NUM_OF_PIXELS_SLICE"];
  materialIDs mateStruct;
  G4double ID;
  //G4int shift_y = 118; //hand
  G4int shift_y = 130; //knee
  //G4int shift_x = 150; //hand
  G4int shift_x = 80; //knee
  for( G4int i_y = 0; i_y < 350 ; i_y++){
	  for( G4int i_x = 0; i_x < 350 ; i_x++){
	    if (i_x < shift_x || i_x > (shift_x + fNVoxelX -1) || i_y < shift_y || i_y > (shift_y + fNVoxelY-1)){
			fin >> ID;
			continue;
	    }
		G4int voxel = (i_x - shift_x) + (i_y - shift_y)*pp->int_map["i_NUM_OF_VOXELS_X"] + voxelCopyNo;

		fin >> ID;
		mateStruct = intensityToMateID[ID];
		mateID = mateStruct.mat1ID;
	//    if (mateID!=0) {
	//    	std::cout << "structureMateID = " << mateID << std::endl;
	//    	std::cout << "ii = " << ii << " voxelCopyNo = " << voxelCopyNo << std::endl;
	//    }
		/*
		//hack to remove non relevant tissue
		//if (i_x > fNVoxelX-15) {mateID =22;};
		*/

		fMateIDs[voxel] = mateID;

//comment when no need to get GT
/*
		output << mateID << "\n";
*/
	  }
  }
//comment when no need to get GT
/*
   output.close();
*/
}



void B1DetectorConstruction::ConstructSDandField()
{
  // -- Fetch volume for biasing:
 // G4LogicalVolume* logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume("phantomContainer");
  //G4LogicalVolume* logicWorld = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  //G4LogicalVolume* logicTestBone = G4LogicalVolumeStore::GetInstance()->GetVolume("bone");

//  // ----------------------------------------------
//  // -- operator creation and attachment to volume:
//  // ----------------------------------------------
//  B1BOptrMultiParticleChangeCrossSection* testMany = new B1BOptrMultiParticleChangeCrossSection();
//  testMany->AddParticle("gamma");
//  testMany->AttachTo(logicTest);
//  testMany->AttachTo(logicWorld);
//   ----------------------------------------------
//   -- operator creation and attachment to volume:
//   ----------------------------------------------
if ((pp->int_map["i_BUILD_PHANTOM"]==1) && (pp->int_map["i_BIASING"]==1) ){
	  B1BOptrFD* FDOptr =  new B1BOptrFD("gamma","FDOperator");
	  FDOptr->AttachTo(fvoxel_logic); //maybe should be attache to fcontainer?
	  //FSOptr->AttachTo(logicWorld);
	  //comptLEOptr->AttachTo(logicTestBone);
	  G4cout << " Attaching biasing operator " << FDOptr->GetName()
			 << " to logical volume " << fvoxel_logic->GetName()
			 << G4endl;
}

	// sensitive detectors
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	//detector1 - SD
	//detector2 - Scorer
	//creating my sensitive detector and adding it to the SD manager - the data will be saved in histograms only if record hist is on
	if (pp->int_map["i_RECORD_HIST"] == 1){
		G4String SDname;
		G4VSensitiveDetector* detector1 = new myDetectorSD(SDname="/detector1");
		SDman->AddNewDetector(detector1);
		//attaching my sensitive detector to the detector logical element
		SetSensitiveDetector(detectorLV,detector1);
	}
	//creating scorer for detector
	G4MultiFunctionalDetector* detector2 = new G4MultiFunctionalDetector("detector2");
	SDman->AddNewDetector(detector2);

	// setting primitive scorers
	for (G4int i=0; i<pp->int_map["i_NUM_OF_SCORERS"]; i++){
		G4VPrimitiveScorer* primitive;
		primitive = new B1EnergyDeposit("eDep_" + IntToString(i),i);
		primitive->SetUnit("keV");
	    detector2->RegisterPrimitive(primitive);
	}

    SetSensitiveDetector(detectorLV,detector2);

}





