/*
 * ParamsSingleton.hh
 *
 *  Created on: Sep 21, 2019
 *      Author: adamgeva
 */

#ifndef PARAMSSINGLETON_HH_
#define PARAMSSINGLETON_HH_

#include <string>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "G4THitsMap.hh"

#include <map>
#include <sstream>
#include "G4Types.hh"
namespace p = boost::python;
namespace np = boost::python::numpy;

class G4run;

class ParamsSingleton{
public:
   static ParamsSingleton* Instance();
   bool parse_dict(p::dict d);
   bool extract_elements_list(p::list list);
   bool extract_A_list(p::list list);
   bool extract_Z_list(p::list list);
   bool extract_materials_array(np::ndarray arr);
   bool extract_dens_array(np::ndarray arr);
   bool extract_ids_array(np::ndarray arr);
   bool extract_src_array(np::ndarray arr);
   bool extract_src_orient_array(np::ndarray arr);
   bool extract_det_array(np::ndarray arr);
   bool extract_det_orient_array(np::ndarray arr);
   bool extract_spectrum_array(np::ndarray arr);
   bool extract_projection_code_array(np::ndarray arr);
   bool extract_intens_code_array(np::ndarray arr);
   bool extract_phantom_offset_orientation_array(np::ndarray arr);
   bool extract_error_array(np::ndarray arr);
   bool create_images_arr();
   bool create_run_images_arr();
   bool create_grad_arr();
   bool create_P_arr();
   bool free_images_mem();
   bool free_run_images_mem();
   bool free_grad_mem();
   bool free_P_mem();
//   bool write_image_run(G4int run, G4int scorer, const G4Run*);
   bool write_grad_run(G4int run, G4double** grad);
   bool write_P_run(G4int run, G4double* P);

   G4double* materials_arr = NULL; // dims: i_NUM_OF_BASE_MATERIALS X i_NUM_OF_ELEMENTS
   G4int* ids_arr = NULL; // dims: N_voxels_Z X N_voxels_Y X N_voxels_X
   G4double* dens_arr = NULL; // dims: N_voxels_Z X N_voxels_Y X N_voxels_X
   G4double* src_arr = NULL; // dims: N_sources X 3
   G4double* src_orient_arr = NULL; // dims: N_sources X 3 X 3
   G4double* det_arr = NULL; // dims: N_rows*N_cols X 3
   G4double* det_orient_arr = NULL; // dims: N_rows*N_cols X 3 X 3
   G4double* spectrum_arr = NULL; // dims: N_spectrum_bins x 1
   G4int* projection_code_arr = NULL; // dims N_shots x N_projections_in_shot
   G4int* intens_code_arr = NULL; // dims N_shots x N_projections_in_shot
   G4double* phantom_loc_orient = NULL; // dims 4 x 3
   G4double* error_arr = NULL; // dims N_shots x N_rows x N_cols
   G4double** run_images_arr = NULL;
   G4double**** images_arr = NULL;
   G4double*** grad_arr = NULL;
   G4double** P_arr = NULL;
   std::vector <std::string> elements;
   std::vector <G4double> A_list;
   std::vector <G4double> Z_list;

   G4int num_phot = 0;


   std::map<std::string, G4double> double_map;
   std::map<std::string, G4int> int_map;
   std::map<std::string, std::string> string_map;


private:
   ParamsSingleton(){};  // Private so that it can  not be called
   void operator=(ParamsSingleton const&){};  // assignment operator is private
   static ParamsSingleton* m_pInstance;


};


#endif /* PARAMSSINGLETON_HH_ */


