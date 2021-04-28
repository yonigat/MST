/*
 * params_dict.cc
 *
 *  Created on: Sep 23, 2019
 *      Author: adamgeva
 */

#include "params_dict.hh"
#include "G4Types.hh"
#include "globalFunctions.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>

p::dict create_params_dict ()
{
	p::dict di;

	di["1"] = 98;
	di["i_MULTI"] = 1;
	di["i_NUM_OF_SCORERS"] = 5;
	di["i_NUM_OF_BASE_SCORERS"] = 5;
	di["i_SINGLE_SCATTERING"] = 0;
	di["i_NUM_OF_SOURCES"] = 1; //this defines the number of runs
	di["i_NUM_OF_SPECTRUM_BINS"] = 150; //from spectrum MATLAB plot
	di["i_NUM_OF_ELEMENTS"] = 6; //this defines the number of elements
	di["i_NUM_OF_THREADS"] = 38;
	di["i_NUM_OF_BASE_MATERIALS"] = 10; // all materials are quantized to 100 materials with different densities
	di["i_NUM_OF_PHOTONS"] = 100000; // equal to mac file
	// verbose
	di["i_VERBOSE_SCORING"] = 0;
	di["i_VERBOSE_PHYSICS_LIST"] = 0;
	di["i_VERBOSE_ANALYSIS_MANAGER"] = 0;
	di["i_VERBOSE_VIS"] = 0;
	di["i_VERBOSE_RUN"] = 1;
	di["i_VERBOSE_EVENT"] = 0;
	di["i_VERBOSE_TRACK"] = 0;
	di["i_VERBOSE_G4NAVIGATOR"] = 0;
	di["i_BIASING"] = 0; //set to 1 for on biasing
	di["i_SPLITTING_FACTOR_NP"] = 70;
	di["i_SPLITTING_FACTOR_NS"] = 70;
	di["f_DETECTOR_SPECIAL_CUTS"] = 1;
	di["f_PHANTOM_PRODUCTION_CUTS"] = 0;
	di["i_KILL_ELECTRONS"] = 1;
	di["i_RECORD_HIST"] = 0; //when 0 - no histograms will be recorded in the simulation
	di["i_PRINT_ELEMENTS_XS"] = 0;
	//Particle Gun
	di["f_PARTICLE_ENERGY"] = 60 ;//KeV
	di["f_MIN_THETA"] = 0;
	di["f_MAX_THETA"] = M_PI/14;
	di["f_MIN_PHI"] = 0;
	di["f_MAX_PHI"] = (2*M_PI);//M_PI/14 ;
	di["f_FILL_FACTOR"] = 0.7;
	di["f_THETA_CUT"] = 6 ;
	di["i_ALT_SOURCES"] = 1;
	//Projection Code
	di["i_NUM_SHOTS"] = 1;
	di["i_NUM_PROJECTION_IN_SHOT"] = 2;
	di["s_FILE_PROJECTION_CODE"] = "../run_inputs/projection_code.txt";
	//Geometry
	di["i_BUILD_DETECTORS"] = 1;
	di["i_CALC_GRADIENT"] = 1;
	di["i_BUILD_PHANTOM"] = 1;
	di["i_XCAT"] = 0;//build XCAT phantom
	di["i_CT_PHANTOM"] = 1; //build phantom using linear reconstruction
	di["s_FILE_XCAT_ID_TO_COMP"] = "../run_inputs/XCAT/id_to_comp.txt";
	di["s_FILE_XCAT_SLICE_PREFIX"] = "../run_inputs/XCAT/hand/out_hand_act_" ;
	di["s_FILE_SOURCES"] = "../run_inputs/src_loc.txt";
	di["s_FILE_SOURCES_ORIENTATION"] = "../run_inputs/src_orient/";
	di["s_FILE_DETECTORS"] = "../run_inputs/det_loc.txt";
	di["s_FILE_DETECTORS_ORIENTATION"] = "../run_inputs/det_orient/";
	di["s_FILE_PHANTOM_OFFSET_ORIENTATION"] = "../run_inputs/phantom_location_orientation.txt";
	di["i_USE_DICOM_INIT"] = 0;
	di["i_ROTATE_SEQUENTIALLY"] = 0;
	di["f_WORLD_XY"] = 210; //cm half size
	di["f_WORLD_Z"] = 210; //cm half size
	di["f_DETECTOR_X"] = 0.8; //mm half size
	di["f_DETECTOR_Y"] = 0.556 ;//mm half size.
	di["f_DETECTOR_Z"] = 0.556; //mm half size. detector depth
	di["i_NUM_OF_DET_COLS"] = 272;
	di["i_NUM_OF_DET_ROWS"] = 224;
	di["f_CENTER_TO_DET"] = 267.38; //mm
	di["f_SOURCE_TO_CENTER"] = 605.27; //mm
	di["f_RADIUS"] = 605.27; //mm
	di["f_OFFSET_U"] = 0; //mm
	di["f_OFFSET_V"] = 0; //mm
	di["f_SHIFT"] = 1; //cm
	di["i_GT_MODE"] = 1;
	di["i_LIV_MODE"] = 1;
	di["f_PHANTOM_OFFSET_X"] = 0; //mm
	di["f_PHANTOM_OFFSET_Y"] = 0; //mm
	di["f_PHANTOM_OFFSET_Z"] = 0; //mm
	//XCAT
	di["i_NUM_OF_Z_SLICES"] = 80;
	di["i_NUM_OF_VOXELS_X"] = 100 ;
	di["i_NUM_OF_VOXELS_Y"] = 100 ;
	di["i_NUM_OF_VOXELS"] = 800000;
	di["i_NUM_OF_PIXELS_SLICE"] = 10000;
	di["f_VOXEL_HALF_X"] = 1.15; //mm
	di["f_VOXEL_HALF_Y"] = 1.15; //mm
	di["f_VOXEL_HALF_Z"] = 1.15; //mm
	di["s_FILE_SPECTRUM"] = "../run_inputs/spectrum.txt";
	di["s_FILE_MATERIALS"] = "../run_inputs/materials.txt"; //v
	di["s_FILE_MATERIALS_DICOM_BASIC"] = "../run_inputs/Dicom_base_materials_new_smart_init.txt";
	di["s_FILE_VOXEL_TO_MATERIALS"] = "../run_inputs/voxels_to_materials.txt";
	di["s_FILE_VOXEL_TO_MATERIALS_ID"] = "../run_inputs/voxels_to_materials_id.txt";
	di["s_FILE_VOXEL_TO_MATERIALS_TEST"] = "../run_inputs/phantom/";
	di["s_FILE_VOXEL_TO_MATERIALS_Y"] = "../run_inputs/Y/";//vx
	di["s_FILE_VOXEL_TO_MATERIALS_DENS"] = "../run_inputs/voxels_to_materials_dens.txt";

	di["s_FILE_FIXED_SOURCE"] = "../run_inputs/all_sources.txt";
	di["s_FILE_SOURCE_TO_DET"] = "../run_outputs_geom/sourceToDet.csv";
	di["s_FILE_SOURCE_POS"] = "../run_outputs_geom/sourcesPos.csv";
	di["s_FILE_DET_POS"] = "../run_outputs_geom/detectorsPos.csv";
	di["s_OUTPUT_DIR"] = "../run_outputs";
	di["s_INPUT_DIR"] = "../run_inputs";
	di["s_GRADIENT_DIR"] = "../run_outputs_grad";
	di["s_ERROR_DIR"] = "../run_inputs/Error/";
	di["s_ELEMENTS_FILE"] = "../run_inputs/elements.txt";
	di["s_A_FILE"] = "../run_inputs/A.txt";
	di["s_Z_FILE"] = "../run_inputs/Z.txt";
	di["i_WRITE_TO_FILES"] = 1;


	return di;
}

p::dict read_params_dict (std::string extracted_val){
	p::dict di;

		std::string fname = extracted_val;
		std::ifstream fin(fname.c_str(), std::ios_base::in);
		if( !fin.is_open() ) {
		   std::cout << "File not found: " + fname << std::endl;
		}
		while (!fin.eof()){
			std::string key;
			fin >> key;
			if (key[0] == 's') {
				std::string str_val;
				fin >> str_val;
				di[key] = str_val;
			} else if (key[0] == 'i') {
				G4int int_val;
				fin >> int_val;
				di[key] = int_val;
			} else if (key[0] == 'f') {
				G4double float_val;
				fin >> float_val;
				di[key] = float_val;
			}

		}

	return di;
}

p::list create_elements_list (p::dict params){
	int numOfEl = p::extract<int> (params["i_NUM_OF_ELEMENTS"]);
	boost::python::list elements_list;

	p::extract<std::string> extracted_val(params["s_ELEMENTS_FILE"]);

	std::string fname = extracted_val;
	std::ifstream fin(fname.c_str(), std::ios_base::in);

	if( !fin.is_open() ) {
	   std::cout << "File not found: " + fname << std::endl;
	}
	for (G4int i=0; i<numOfEl; i++){
		if( fin.eof() ) break;
		std::string el;
		fin >> el;
		elements_list.append(el);

	}
	return elements_list;
}

p::list create_Z_list (p::dict params){
	int numOfEl = p::extract<int> (params["i_NUM_OF_ELEMENTS"]);
	boost::python::list Z_list;

	p::extract<std::string> extracted_val(params["s_Z_FILE"]);
	std::string fname = extracted_val;
	std::ifstream fin(fname.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + fname << std::endl;
	}
	for (G4int i=0; i<numOfEl; i++){
		if( fin.eof() ) break;
		G4double z;
		fin >> z;
		Z_list.append(z);

	}
	return Z_list;
}

p::list create_A_list (p::dict params){
	int numOfEl = p::extract<int> (params["i_NUM_OF_ELEMENTS"]);
	boost::python::list A_list;

	p::extract<std::string> extracted_val(params["s_A_FILE"]);
	std::string fname = extracted_val;
	std::ifstream fin(fname.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + fname << std::endl;
	}
	for (G4int i=0; i<numOfEl; i++){
		if( fin.eof() ) break;
		G4double a;
		fin >> a;
		A_list.append(a);

	}
	return A_list;
}


np::ndarray create_material_array (p::dict params){
	// Cretae the array:
	int numOfEl = p::extract<int> (params["i_NUM_OF_ELEMENTS"]);

	int numOfBaseMaterials = p::extract<int> (params["i_NUM_OF_BASE_MATERIALS"]);

	p::tuple shape = p::make_tuple(numOfBaseMaterials, numOfEl);
	np::dtype dtype = np::dtype::get_builtin<G4double>();
	np::ndarray a = np::zeros(shape, dtype);

	//Read from file:
	p::extract<std::string> extracted_val(params["s_FILE_MATERIALS"]);
	std::string fname = extracted_val;
	std::ifstream fin(fname.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + fname << std::endl;
	}
	G4int mat;
	G4double cell;
	for (mat = 0; mat<numOfBaseMaterials; mat ++){
		if( fin.eof() ) break;
		// read fractions
		for (G4int i=0; i<numOfEl; i++){
			fin >> cell;
			a[mat][i] = cell;
		}
	}

	std::cout << "Original array:\n" << p::extract<char const *>(p::str(a)) << std::endl;
	fin.close();
	return a;
}

np::ndarray create_ids_array (p::dict params){

	int fNVoxelX = p::extract<int>(params["i_NUM_OF_VOXELS_X"]);
	int fNVoxelY = p::extract<int>(params["i_NUM_OF_VOXELS_Y"]);
	int fNVoxelZ = p::extract<int>(params["i_NUM_OF_Z_SLICES"]);
	p::extract<std::string> extracted_val(params["s_FILE_VOXEL_TO_MATERIALS_Y"]);

	p::tuple shape = p::make_tuple(fNVoxelZ, fNVoxelY, fNVoxelX);
	np::dtype dtype = np::dtype::get_builtin<G4int>();
	np::ndarray ids_arr = np::zeros(shape, dtype);

	std::string file = extracted_val;
	for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
		std::string fileName_id = file + "id/id" + IntToString(iz) + ".txt";
		std::ifstream fin_id(fileName_id.c_str(), std::ios_base::in);
		if( !fin_id.is_open() ) {
		   std::cout << "File not found: " + fileName_id << std::endl;
		}

//		std::cout << "reading file: " << fileName_id << std::endl;
		for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
			  for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
				G4int mateID;
				fin_id >> mateID;
				ids_arr[iz][iy][ix] = mateID;
			  }
		}
		fin_id.close();
	}
	return ids_arr;

}

np::ndarray create_dens_array (p::dict params){

	int fNVoxelX = p::extract<int>(params["i_NUM_OF_VOXELS_X"]);
	int fNVoxelY = p::extract<int>(params["i_NUM_OF_VOXELS_Y"]);
	int fNVoxelZ = p::extract<int>(params["i_NUM_OF_Z_SLICES"]);
	p::extract<std::string> extracted_val(params["s_FILE_VOXEL_TO_MATERIALS_Y"]);

	p::tuple shape = p::make_tuple(fNVoxelZ, fNVoxelY, fNVoxelX);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray dens_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	for( G4int iz = 0; iz < fNVoxelZ; iz++ ) {
		std::string fileName_dens = file + "dens/dens" + IntToString(iz) + ".txt";
		std::ifstream fin_dens(fileName_dens.c_str(), std::ios_base::in);
		if( !fin_dens.is_open() ) {
		   std::cout << "File not found: " + fileName_dens << std::endl;
		}
//		std::cout << "reading file: " << fileName_dens << std::endl;
		for( G4int iy = 0; iy < fNVoxelY; iy++ ) {
			  for( G4int ix = 0; ix < fNVoxelX; ix++ ) {
				G4double dens;
				fin_dens >> dens;
				dens_arr[iz][iy][ix] = dens;
			  }
		}
		fin_dens.close();
	}
	return dens_arr;

}

np::ndarray create_src_array (p::dict params){
	int iNSrcs = p::extract<int>(params["i_NUM_OF_SOURCES"]);
	int iDims = 3;
	p::extract<std::string> extracted_val(params["s_FILE_SOURCES"]);

	p::tuple shape = p::make_tuple(iNSrcs, iDims);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray det_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	std::ifstream fin(file.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}
	G4int src;
	G4double cell;
		for (src = 0; src<iNSrcs; src ++){
			if( fin.eof() ) break;
			// read locations
			for (G4int dim=0; dim<iDims; dim++){
				fin >> cell;
				det_arr[src][dim] = cell;
			}
		}

	return det_arr;
}

np::ndarray create_src_orient_array (p::dict params){

	int iNSrcs = p::extract<int>(params["i_NUM_OF_SOURCES"]);
	int iDims = 3;
	p::extract<std::string> extracted_val(params["s_FILE_SOURCES_ORIENTATION"]);

	p::tuple shape = p::make_tuple(iNSrcs, iDims, iDims);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray src_orient_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
//	std::string fileName_dens = file + "src_orient" + IntToString(src) + ".txt";
	std::ifstream fin(file.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}
	G4double cell;
	for( G4int src = 0; src < iNSrcs; src++ ) {
		for( G4int iy = 0; iy < iDims; iy++ ) {
			  for( G4int ix = 0; ix < iDims; ix++ ) {
				fin >> cell;
				src_orient_arr[src][iy][ix] = cell;
			  }
		}
	}
	fin.close();
	return src_orient_arr;

}
np::ndarray create_det_array (p::dict params){
	int iNRows = p::extract<int>(params["i_NUM_OF_DET_ROWS"]);
	int iNCols = p::extract<int>(params["i_NUM_OF_DET_COLS"]);
	int iDets = iNRows * iNCols;
	int iDims = 3;
	p::extract<std::string> extracted_val(params["s_FILE_DETECTORS"]);

	p::tuple shape = p::make_tuple(iDets, iDims);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray det_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	std::ifstream fin(file.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}
	G4int det;
	G4double cell;
		for (det = 0; det<iDets; det ++){
			if( fin.eof() ) break;
			// read locations
			for (G4int dim=0; dim<iDims; dim++){
				fin >> cell;
				det_arr[det][dim] = cell;
			}
		}

	return det_arr;
}

np::ndarray create_det_orient_array (p::dict params){

	int iNRows = p::extract<int>(params["i_NUM_OF_DET_ROWS"]);
	int iNCols = p::extract<int>(params["i_NUM_OF_DET_COLS"]);
	int iDets = iNRows * iNCols;
	int iDims = 3;
	p::extract<std::string> extracted_val(params["s_FILE_DETECTORS_ORIENTATION"]);

	p::tuple shape = p::make_tuple(iDets, iDims, iDims);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray det_orient_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;

//	std::string fileName_dens = file + "det_orient" + IntToString(det) + ".txt";
	std::ifstream fin_dens(file.c_str(), std::ios_base::in);
	if( !fin_dens.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}
	G4double dens;
	for( G4int det = 0; det < iDets; det++ ) {
		for( G4int iy = 0; iy < iDims; iy++ ) {
			  for( G4int ix = 0; ix < iDims; ix++ ) {
				fin_dens >> dens;
				det_orient_arr[det][iy][ix] = dens;
			  }
		}
	}
	fin_dens.close();
	return det_orient_arr;

}

np::ndarray create_spectrum_array (p::dict params){
	int iNBins = p::extract<int>(params["i_NUM_OF_SPECTRUM_BINS"]);

	p::extract<std::string> extracted_val(params["s_FILE_SPECTRUM"]);

	p::tuple shape = p::make_tuple(iNBins);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray spectrum_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	std::ifstream fin(file.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}

	G4double cell;
		for (G4int bin = 0; bin<iNBins; bin ++){
			if( fin.eof() ) break;
			// read locations
			fin >> cell;
			spectrum_arr[bin] = cell;
		}

	return spectrum_arr;
}

np::ndarray create_projection_code_array (p::dict params){
	int iNShots = p::extract<int>(params["i_NUM_SHOTS"]);
	int iNProject = p::extract<int>(params["i_NUM_PROJECTION_IN_SHOT"]);

	p::extract<std::string> extracted_val(params["s_FILE_PROJECTION_CODE"]);

	p::tuple shape = p::make_tuple(iNShots, iNProject);
	np::dtype dtype_double = np::dtype::get_builtin<G4int>();
	np::ndarray code_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	std::ifstream fin(file.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}

	G4int cell;
		for (G4int shot = 0; shot<iNShots; shot ++){
			if( fin.eof() ) break;
			// read locations
			for (G4int project=0; project<iNProject; project++){
				fin >> cell;
				code_arr[shot][project] = cell;
			}
		}

	return code_arr;
}

np::ndarray create_intens_code_array (p::dict params){
	int iNShots = p::extract<int>(params["i_NUM_SHOTS"]);
	int iNProject = p::extract<int>(params["i_NUM_PROJECTION_IN_SHOT"]);

	p::extract<std::string> extracted_val(params["s_FILE_INTENS_CODE"]);

	p::tuple shape = p::make_tuple(iNShots, iNProject);
	np::dtype dtype_double = np::dtype::get_builtin<G4int>();
	np::ndarray code_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	std::ifstream fin(file.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}

	G4int cell;
		for (G4int shot = 0; shot<iNShots; shot ++){
			if( fin.eof() ) break;
			// read locations
			for (G4int project=0; project<iNProject; project++){
				fin >> cell;
				code_arr[shot][project] = cell;
			}
		}

	return code_arr;
}

np::ndarray create_phantom_loc_orient_array (p::dict params){
	int iNCols = 3;
	int iNRows = 4;
	p::extract<std::string> extracted_val(params["s_FILE_PHANTOM_OFFSET_ORIENTATION"]);

	p::tuple shape = p::make_tuple(iNRows, iNCols);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray det_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	std::ifstream fin(file.c_str(), std::ios_base::in);
	if( !fin.is_open() ) {
	   std::cout << "File not found: " + file << std::endl;
	}
	G4int row;
	G4double cell;
		for (row = 0; row<iNRows; row ++){
			if( fin.eof() ) break;
			// read locations
			for (G4int col=0; col<iNCols; col++){
				fin >> cell;
				det_arr[row][col] = cell;
			}
		}

	return det_arr;
}

np::ndarray create_error_array (p::dict params){

	int iNRuns = p::extract<int>(params["i_NUM_SHOTS"]);
	int iNRows = p::extract<int>(params["i_NUM_OF_DET_ROWS"]);
	int iNCols = p::extract<int>(params["i_NUM_OF_DET_COLS"]);
	p::extract<std::string> extracted_val(params["s_ERROR_DIR"]);

	p::tuple shape = p::make_tuple(iNRuns, iNCols, iNRows);
	np::dtype dtype_double = np::dtype::get_builtin<G4double>();
	np::ndarray error_arr = np::zeros(shape, dtype_double);

	std::string file = extracted_val;
	for( G4int run = 0; run < iNRuns; run++ ) {
		std::string fileName_error = file + "error_" + IntToString(run) + ".txt";
		std::ifstream fin_dens(fileName_error.c_str(), std::ios_base::in);
		if( !fin_dens.is_open() ) {
		   std::cout << "File not found: " + fileName_error << std::endl;
		}
//		std::cout << "reading file: " << fileName_error << std::endl;
		for( G4int iy = 0; iy < iNCols; iy++ ) {
			  for( G4int ix = 0; ix < iNRows; ix++ ) {
				G4double error;
				fin_dens >> error;
				error_arr[run][iy][ix] = error;
			  }
		}
		fin_dens.close();
	}
	return error_arr;

}


