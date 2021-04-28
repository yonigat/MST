/*
 * params_dict.hh
 *
 *  Created on: Sep 23, 2019
 *      Author: adamgeva
 */

#ifndef PARAMS_DICT_HH_
#define PARAMS_DICT_HH_

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace np = boost::python::numpy;


p::dict create_params_dict ();
p::dict read_params_dict (std::string extracted_val);
p::list create_elements_list (p::dict params);
p::list create_A_list (p::dict params);
p::list create_Z_list (p::dict params);
np::ndarray create_material_array (p::dict params);
np::ndarray create_ids_array (p::dict params);
np::ndarray create_dens_array (p::dict params);
np::ndarray create_src_array (p::dict params);
np::ndarray create_src_orient_array (p::dict params);
np::ndarray create_det_array (p::dict params);
np::ndarray create_det_orient_array (p::dict params);
np::ndarray create_spectrum_array (p::dict params);
np::ndarray create_projection_code_array (p::dict params);
np::ndarray create_intens_code_array (p::dict params);
np::ndarray create_phantom_loc_orient_array (p::dict params);
np::ndarray create_error_array (p::dict params);


#endif /* PARAMS_DICT_HH_ */
