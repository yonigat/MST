/*
 * ParamsSingleton.cc
 *
 *  Created on: Sep 21, 2019
 *      Author: adamgeva
 */
#include "ParamsSingleton.hh"
#include <stddef.h>  // defines NULL
#include <iostream>
#include "G4THitsMap.hh"
#include "B1Run.hh"
namespace p = boost::python;
namespace np = boost::python::numpy;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
*/

// Global static pointer used to ensure a single instance of the class.
ParamsSingleton* ParamsSingleton::m_pInstance = NULL;

ParamsSingleton* ParamsSingleton::Instance()
{
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		  m_pInstance = new ParamsSingleton;
//		  m_pInstance->images_arr = new G4double;
	return m_pInstance;
}

bool ParamsSingleton::parse_dict(boost::python::dict d)
{

	std::string s_value;
	float f_value;
	int i_value;

	boost::python::list keys = d.keys();
	for (int i = 0; i < len(keys); ++i) {
	   boost::python::extract<std::string> extracted_key(keys[i]);
	   if(!extracted_key.check()){
			std::cout<<"Key invalid, map might be incomplete"<<std::endl;
			continue;
	   }
	   std::string key = extracted_key;
	   if (key[0]=='i'){
		   boost::python::extract<G4int> extracted_val(d[key]);
		   if(!extracted_val.check()){
		   std::cout<<"Value invalid, map might be incomplete"<<std::endl;
				continue;
		   }
		   i_value = extracted_val;
		   int_map[key] = i_value;
	   } else if (key[0]=='f'){
		   boost::python::extract<G4double> extracted_val(d[key]);
		   if(!extracted_val.check()){
		   std::cout<<"Value invalid, map might be incomplete"<<std::endl;
				continue;
		   }
		   f_value = extracted_val;
		   double_map[key] = f_value;
	   } else if (key[0]=='s'){
		   boost::python::extract<std::string> extracted_val(d[key]);
		   if(!extracted_val.check()){
		   std::cout<<"Value invalid, map might be incomplete"<<std::endl;
				continue;
		   }
		   s_value = extracted_val;
		   string_map[key] = s_value;
	   }
	 }
	 return true;
}

bool ParamsSingleton::extract_elements_list(p::list list)
{
	for (int i = 0; i < len(list); ++i) {
		   boost::python::extract<std::string> extracted_el(list[i]);
		   std::string el = extracted_el;
		   elements.push_back(el);
	}
	return true;
}

bool ParamsSingleton::extract_Z_list(p::list list)
{
	for (int i = 0; i < len(list); ++i) {
		   boost::python::extract<G4double> extracted_val(list[i]);
		   G4double z = extracted_val;
		   Z_list.push_back(z);
//		   std::cout << "element " << i << " " << Z_list[i] << std::endl;
	}
	return true;
}

bool ParamsSingleton::extract_A_list(p::list list)
{
	for (int i = 0; i < len(list); ++i) {
		   boost::python::extract<G4double> extracted_val(list[i]);
		   G4double a = extracted_val;
		   A_list.push_back(a);
//		   std::cout << "element " << i << " " << A_list[i] << std::endl;
	}
	return true;
}

bool ParamsSingleton::extract_materials_array(np::ndarray arr)
{
	materials_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_dens_array(np::ndarray arr)
{
	dens_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_ids_array(np::ndarray arr)
{
	ids_arr = reinterpret_cast <G4int *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_src_array(np::ndarray arr)
{
	src_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_src_orient_array(np::ndarray arr)
{
	src_orient_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_det_array(np::ndarray arr)
{
	det_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_det_orient_array(np::ndarray arr)
{
	det_orient_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_spectrum_array(np::ndarray arr)
{
	spectrum_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_projection_code_array(np::ndarray arr)
{
	projection_code_arr = reinterpret_cast <G4int *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_intens_code_array(np::ndarray arr)
{
	intens_code_arr = reinterpret_cast <G4int *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_phantom_offset_orientation_array(np::ndarray arr)
{
	phantom_loc_orient = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::extract_error_array(np::ndarray arr)
{
	error_arr = reinterpret_cast <G4double *> (arr.get_data());
	return true;
}

bool ParamsSingleton::create_grad_arr()
{
	grad_arr = new G4double**[int_map["i_NUM_SHOTS"]];
	for (G4int i=0; i<int_map["i_NUM_SHOTS"]; i++){
		grad_arr[i] = new G4double*[int_map["i_NUM_OF_VOXELS"]];
		for (G4int j=0; j<int_map["i_NUM_OF_VOXELS"]; j++){
			grad_arr[i][j] = new G4double[int_map["i_NUM_OF_ELEMENTS"]];
		}
	}
	return true;
}

bool ParamsSingleton::create_P_arr()
{
	P_arr = new G4double*[int_map["i_NUM_SHOTS"]];
	for (G4int i=0; i<int_map["i_NUM_SHOTS"]; i++){
		P_arr[i] = new G4double[int_map["i_NUM_OF_VOXELS"]];
	}
	return true;
}

bool ParamsSingleton::create_run_images_arr()
{
	run_images_arr = new G4double*[int_map["i_NUM_OF_SCORERS"]];
	for (G4int i=0; i<int_map["i_NUM_OF_SCORERS"]; i++){
		run_images_arr[i] = new G4double[int_map["i_NUM_OF_DET_COLS"]*int_map["i_NUM_OF_DET_ROWS"]];
	}
	return true;
}

bool ParamsSingleton::create_images_arr()
{
	images_arr = new G4double***[int_map["i_NUM_SHOTS"]];
	for (G4int i=0; i<int_map["i_NUM_SHOTS"]; i++){
		images_arr[i] = new G4double**[int_map["i_NUM_OF_SCORERS"]];
		for (G4int j=0; j<int_map["i_NUM_OF_SCORERS"]; j++){
			images_arr[i][j] = new G4double*[int_map["i_NUM_OF_DET_ROWS"]];
			for (G4int k=0; k<int_map["i_NUM_OF_DET_ROWS"]; k++){
				images_arr[i][j][k] = new G4double[int_map["i_NUM_OF_DET_COLS"]];
			}
		}
	}
//	for (G4int i=0; i<int_map["i_NUM_SHOTS"]; i++){
//		for (G4int j=0; j<int_map["i_NUM_OF_SCORERS"]; j++){
//			for (G4int k=0; k<int_map["i_NUM_OF_DET_ROWS"]; k++){
//				for (G4int l=0; l<int_map["i_NUM_OF_DET_COLS"]; l++){
//				images_arr[i][j][k][l] = 0;
//				}
//			}
//		}
//	}
	return true;
}

bool ParamsSingleton::free_images_mem(){
	for (G4int i=0; i<int_map["i_NUM_SHOTS"]; i++){
		for (G4int j=0; j<int_map["i_NUM_OF_SCORERS"]; j++){
			for (G4int k=0; k<int_map["i_NUM_OF_DET_ROWS"]; k++)
				delete [] images_arr[i][j][k];
			delete [] images_arr[i][j];
		}
		delete [] images_arr[i];
	}
	delete [] images_arr;
	return true;
}

bool ParamsSingleton::free_P_mem(){
	for (G4int i=0; i<int_map["i_NUM_SHOTS"]; i++){
		delete [] P_arr[i];
	}
	delete [] P_arr;
	return true;
}

bool ParamsSingleton::free_run_images_mem(){
	for (G4int i=0; i<int_map["i_NUM_OF_SCORERS"]; i++){
		delete [] run_images_arr[i];
	}
	delete [] run_images_arr;
	return true;
}

bool ParamsSingleton::free_grad_mem(){
	for (G4int i=0; i<int_map["i_NUM_SHOTS"]; i++){
		for (G4int j=0; j<int_map["i_NUM_OF_VOXELS"]; j++)
			delete [] grad_arr[i][j];
		delete [] grad_arr[i];
	}
	delete [] grad_arr;
	return true;
}

bool ParamsSingleton::write_grad_run(G4int run, G4double** grad){
	for (G4int i=0; i<int_map["i_NUM_OF_VOXELS"]; i++){
		for (G4int j=0; j<int_map["i_NUM_OF_ELEMENTS"]; j++){
			grad_arr[run][i][j] = grad[i][j];
		}
	}
	return true;
}

bool ParamsSingleton::write_P_run(G4int run, G4double* P){
	for (G4int i=0; i<int_map["i_NUM_OF_VOXELS"]; i++){
			P_arr[run][i] = P[i];
	}
	return true;
}

/*
bool ParamsSingleton::write_image_run(G4int run, G4int scorer, const G4Run* aRun){
	const B1Run* theRun = (const B1Run*)aRun;
	for (G4int i=0; i<int_map["i_NUM_OF_ROWS"]; i++){
		for (G4int j=0; j<int_map["i_NUM_OF_COLS"]; j++){
			G4double* edep =  (theRun->fMapSum[scorer])[i*int_map["i_NUM_OF_DET_ROWS"] +j];
			if (edep == 0)
				images_arr[run][scorer][i][j] = 0;
			else
				images_arr[run][scorer][i][j] = *edep;
		}
	}
	return true;
}
*/
