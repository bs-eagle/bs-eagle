/** 
 * @file py_dens.h
 * @brief Python interface for dens_matrix 
 * @date 2009-12-09
 */
#ifndef __PY_DENS_H
#define __PY_DENS_H

#include "dens_matrix_iface.h"
#include "dens_matrix.h"
#include "dens_matrix_tools.h"
#include "py_matrix_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

    // Sergey Miryanov at 07.04.2008
    // Refactored at 16.10.2009
    // for casting into python from child classes to parent class
    // we should export base class
#if 0    
  /** 
   * @brief python wrapper for dens_matrix_iface class
   * 
   */
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), dens_matrix_iface, py_dens_matrix_iface)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <fp_vector_type, i_vector_type>        matrix_t;               
    typedef dens_matrix_iface <fp_vector_type, i_vector_type, fp_storage_vector_type>    dens_matrix_t;         
                                                                                        
    typedef typename fp_vector_type::value_type                 fp_type_t;              
    typedef typename fp_storage_vector_type::value_type         fp_storage_type_t;      
    typedef typename i_vector_type::value_type                  i_type_t;               
                                                                                        
                                                                                        
    //typedef dens_matrix_iface <fp_vector_type, i_vector_type, fp_storage_vector_type>         wrapped_t;              
    //typedef BOOST_PP_CAT (py_dens_matrix_iface, _base) <fp_vector_type, i_vector_type, fp_storage_vector_type>        base_t;                 
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), py_dens_matrix_iface);
  
  public:
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (calc_lin_comb, int, 5, (fp_type_t, fp_type_t, const fp_vector_type&, const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (get_allocated_memory_in_mbytes, fp_type_t, 0, (empty_arg__));       
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_block_size, i_type_t, 0, (empty_arg__));            
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_rows, i_type_t, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_cols, i_type_t, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_CONST   (init_vector, void, 1, (fp_vector_type &));                       
    WRAP_PURE_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                   
    WRAP_PURE_METHOD_R       (init_by_matrix, int, 1, (const dens_matrix_t&));          
    WRAP_PURE_METHOD_R       (init, int, 3, (const i_type_t, const i_type_t, const i_type_t));          
    WRAP_PURE_METHOD_R       (copy, int, 1, (const dens_matrix_t&));                    
    WRAP_PURE_METHOD_R_CONST (get_calc_block_size, i_type_t, 0, (empty_arg__));             
    WRAP_PURE_METHOD         (set_calc_block_size, void, 1, (const i_type_t));             
    WRAP_PURE_METHOD_R       (get_values, fp_storage_vector_type&, 0, (empty_arg__));   
    WRAP_PURE_METHOD_R_CONST (get_values_const, const fp_storage_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R_CONST (internal_check, int, 0, (empty_arg__));                   
  };
#endif //0

  PY_EXPORTER (py_dens_matrix_iface_exporter, py_matrix_iface_exporter)
    .add_property ("calc_block_size",           
                   &T::get_calc_block_size, &T::set_calc_block_size, 
                   "Block size for calculation algorithms")
    .def ("init_by_matrix",                     
          &T::init_by_matrix,
          args ("matrix"), "Initialize matrix by matrix :)")
    .def ("init",                               
          &T::init,
          args ("n_rows", "n_cols", "calc_block_size"), "Initialize matrix")
    .def ("copy",                               
          &T::copy,
          args ("matrix"), "Initialize matrix by matrix and copy all content")
    .def ("get_values",                         
          &T::get_values, //return_value_policy <reference_existing_object> (),
          args (""), "Return reference to values vector")
    //.def ("get_values_const",                   
    //      &T::get_values_const, return_value_policy <reference_existing_object> (),
    //      args (""), "Return const reference to values vector")
    .def ("internal_check",                     
          &T::internal_check,
          args (""), "Return 0 if OK")
  PY_EXPORTER_END;                               

  PY_EXPORTER (py_matrix_dens_tools_exporter, default_exporter)
    .def ("random_init",                        
          &T::random_init,
          args ("matrix", "n_rows", "n_cols", "calc_block_size", "value_dispertion"),
          "Initialize matrix by random values")
  PY_EXPORTER_END;                               

  //! export matrices to python                  
  void py_export_dens_matrices ();                    

    
  } // namespace python
} // namespace blue_sky
#endif //

#endif //__PY_DENS_H

