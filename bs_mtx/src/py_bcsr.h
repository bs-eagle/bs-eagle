/** 
 * @file py_bcsr.h
 * @brief python interface for BCSR matrix
 * @author 
 * @date 2009-12-07
 */
#ifndef __PY_BCSR_H
#define __PY_BCSR_H

#include "bcsr_matrix_iface.h"
#include "bcsr_amg_matrix_iface.h"
#include "bcsr.h"
#include "bcsr_matrix_tools.h"
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
 * @brief python wrapper for bcsr_matrix_iface class
 * 
 */
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), bcsr_matrix_iface, py_bcsr_matrix_iface)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <fp_vector_type, i_vector_type>        matrix_t;               
    typedef bcsr_matrix_iface <fp_vector_type, i_vector_type, fp_storage_vector_type>    bcsr_matrix_t;         
                                                                                        
    typedef typename fp_vector_type::value_type                 t_double;              
    typedef typename fp_storage_vector_type::value_type         t_float;      
    typedef typename i_vector_type::value_type                  t_long;               
                                                                                        
                                                                                        
    typedef bcsr_matrix_iface <fp_vector_type, i_vector_type, fp_storage_vector_type>         wrapped_t;              
    typedef BOOST_PP_CAT (py_bcsr_matrix_iface, _base) <fp_vector_type, i_vector_type, fp_storage_vector_type>        base_t;                 
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), py_bcsr_matrix_iface);
  
  public:
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (calc_lin_comb, int, 5, (t_double, t_double, const fp_vector_type&, const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (get_allocated_memory_in_mbytes, t_double, 0, (empty_arg__));       
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_block_size, t_long, 0, (empty_arg__));            
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_rows, t_long, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_cols, t_long, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_CONST   (init_vector, void, 1, (fp_vector_type &));                       
    WRAP_PURE_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                   
    WRAP_PURE_METHOD_R       (init_by_matrix, int, 1, (const bcsr_matrix_t&));          
    WRAP_PURE_METHOD_R       (init, int, 4, (const t_long, const t_long, const t_long, const t_long));          
    WRAP_PURE_METHOD_R       (init_struct, int, 3, (const t_long, const t_long, const t_long));          
    WRAP_PURE_METHOD_R       (alloc_rows_ptr, int, 2, (const t_long, t_long));      
    WRAP_PURE_METHOD_R       (alloc_cols_ind, int, 1, (const t_long));                
    WRAP_PURE_METHOD_R       (alloc_values, int, 2, (const t_long, const t_long));  
    WRAP_PURE_METHOD_R       (alloc_cols_ind_and_values, int, 2, (const t_long, const t_long));  
    WRAP_PURE_METHOD         (set_n_cols, void, 1, (const t_long));                   
    WRAP_PURE_METHOD_R       (copy, int, 1, (const bcsr_matrix_t&));                    
    WRAP_PURE_METHOD_R_CONST (get_n_non_zeros, t_long, 0, (empty_arg__));             
    WRAP_PURE_METHOD_R       (get_rows_ptr, i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R_CONST (get_rows_ptr_const, const i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R       (get_cols_ind, i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R_CONST (get_cols_ind_const, const i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R       (get_values, fp_storage_vector_type&, 0, (empty_arg__));   
    WRAP_PURE_METHOD_R_CONST (get_values_const, const fp_storage_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R_CONST (internal_check, int, 0, (empty_arg__));                   
  };

  /** 
   * @brief python wrapper for bcsr_amg_matrix_iface class
   * 
   */
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), bcsr_amg_matrix_iface, py_bcsr_amg_matrix_iface)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <fp_vector_type, i_vector_type>        matrix_t;               
    typedef bcsr_matrix_iface <fp_vector_type, i_vector_type, fp_storage_vector_type>    bcsr_matrix_t;         
                                                                                        
    typedef typename fp_vector_type::value_type                 t_double;              
    typedef typename fp_storage_vector_type::value_type         t_float;      
    typedef typename i_vector_type::value_type                  t_long;               
                                                                                        
                                                                                        
    typedef bcsr_amg_matrix_iface <fp_vector_type, i_vector_type, fp_storage_vector_type>         wrapped_t;              
    typedef BOOST_PP_CAT (py_bcsr_amg_matrix_iface, _base) <fp_vector_type, i_vector_type, fp_storage_vector_type>        base_t;                 
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), py_bcsr_amg_matrix_iface);
  
  public:
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (calc_lin_comb, int, 5, (t_double, t_double, const fp_vector_type&, const fp_vector_type&, fp_vector_type&)); 
    WRAP_PURE_METHOD_R_CONST (get_allocated_memory_in_mbytes, t_double, 0, (empty_arg__));       
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_block_size, t_long, 0, (empty_arg__));            
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_rows, t_long, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_cols, t_long, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_CONST   (init_vector, void, 1, (fp_vector_type &));                       
    WRAP_PURE_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                   
    WRAP_PURE_METHOD_R       (init_by_matrix, int, 1, (const bcsr_matrix_t&));          
    WRAP_PURE_METHOD_R       (init, int, 4, (const t_long, const t_long, const t_long, const t_long));          
    WRAP_PURE_METHOD_R       (init_struct, int, 3, (const t_long, const t_long, const t_long));          
    WRAP_PURE_METHOD_R       (alloc_rows_ptr, int, 2, (const t_long, t_long));      
    WRAP_PURE_METHOD_R       (alloc_cols_ind, int, 1, (const t_long));                
    WRAP_PURE_METHOD_R       (alloc_values, int, 2, (const t_long, const t_long));  
    WRAP_PURE_METHOD_R       (alloc_cols_ind_and_values, int, 2, (const t_long, const t_long));  
    WRAP_PURE_METHOD         (set_n_cols, void, 1, (const t_long));                    
    WRAP_PURE_METHOD_R       (copy, int, 1, (const bcsr_matrix_t&));                    
    WRAP_PURE_METHOD_R       (build_transpose, int, 4, (const bcsr_matrix_t&, const t_long, const t_long, const t_long));  
    WRAP_PURE_METHOD_R       (build_transpose_struct, int, 4, (const bcsr_matrix_t&, const t_long, const t_long, const t_long));  
    WRAP_PURE_METHOD_R       (build_transpose_struct, int, 4, (const bcsr_matrix_t&, const bcsr_matrix_t&, const bcsr_matrix_t&, const t_long));  
    WRAP_PURE_METHOD_R_CONST (get_n_non_zeros, t_long, 0, (empty_arg__));             
    WRAP_PURE_METHOD_R       (get_rows_ptr, i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R_CONST (get_rows_ptr_const, const i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R       (get_cols_ind, i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R_CONST (get_cols_ind_const, const i_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R       (get_values, fp_storage_vector_type&, 0, (empty_arg__));   
    WRAP_PURE_METHOD_R_CONST (get_values_const, const fp_storage_vector_type&, 0, (empty_arg__));          
    WRAP_PURE_METHOD_R_CONST (internal_check, int, 0, (empty_arg__));                   
  };

      /** 
       * @brief python wrapper for bcsr_matrix_iface class
       */
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), bcsr, py_bcsr)
  {                                                                                             
  public:                                                                                       
    typedef matrix_iface <fp_vector_type, i_vector_type>        matrix_t;                       
    typedef bcsr_matrix_iface <fp_vector_type, i_vector_type, fp_storage_vector_type>    bcsr_matrix_t;         
                                                                                                
    typedef typename fp_vector_type::value_type                 t_double;                      
    typedef typename fp_storage_vector_type::value_type         t_float;              
    typedef typename i_vector_type::value_type                  t_long;                       
                                                                                                
                                                                                                
    typedef bcsr <fp_vector_type, i_vector_type, fp_storage_vector_type>         wrapped_t;    
    typedef BOOST_PP_CAT (py_bcsr, _base) <fp_vector_type, i_vector_type, fp_storage_vector_type>        base_t;                 
                                                                                                
  public:                                                                                       
    CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), py_bcsr);

  public:
    WRAPPER_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAPPER_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type&, fp_vector_type&)); 
    WRAPPER_METHOD_R_CONST (calc_lin_comb, int, 5, (t_double, t_double, const fp_vector_type&, const fp_vector_type&, fp_vector_type&)); 
    WRAPPER_METHOD_R_CONST (get_allocated_memory_in_mbytes, t_double, 0, (empty_arg__));         
    WRAPPER_METHOD_R_CONST (get_n_block_size, t_long, 0, (empty_arg__));                      
    WRAPPER_METHOD_R_CONST (get_n_rows, t_long, 0, (empty_arg__));                            
    WRAPPER_METHOD_R_CONST (get_n_cols, t_long, 0, (empty_arg__));                            
    WRAPPER_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                                 
    WRAPPER_METHOD_CONST   (init_vector, void, 1, (fp_vector_type &));                       
    WRAPPER_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                             
    WRAPPER_METHOD_R       (init_by_matrix, int, 1, (const bcsr_matrix_t&));                    
    WRAPPER_METHOD_R       (init, int, 4, (const t_long, const t_long, const t_long, const t_long));          
    WRAPPER_METHOD_R       (init_struct, int, 3, (const t_long, const t_long, const t_long));          
    WRAPPER_METHOD_R       (alloc_rows_ptr, int, 2, (const t_long, t_long));                
    WRAPPER_METHOD_R       (alloc_cols_ind, int, 1, (const t_long));                          
    WRAPPER_METHOD_R       (alloc_values, int, 2, (const t_long, const t_long));            
    WRAPPER_METHOD_R       (alloc_cols_ind_and_values, int, 2, (const t_long, const t_long));  
    WRAPPER_METHOD         (set_n_cols, void, 1, (const t_long));                              
    WRAPPER_METHOD_R       (init_diag_ind, int, 2, (int, t_long));                            
    WRAPPER_METHOD_R       (copy, int, 1, (const bcsr_matrix_t&));                              
    WRAPPER_METHOD_R       (build_transpose, int, 4, (const bcsr_matrix_t&, const t_long, const t_long, const t_long));  
    WRAPPER_METHOD_R       (build_transpose_struct, int, 4, (const bcsr_matrix_t&, const t_long, const t_long, const t_long));  
    WRAPPER_METHOD_R       (build_transpose_struct, int, 4, (const bcsr_matrix_t&, const bcsr_matrix_t&, const bcsr_matrix_t&, const t_long));  
    WRAPPER_METHOD_R_CONST (get_n_non_zeros, t_long, 0, (empty_arg__));                       
    WRAPPER_METHOD_R       (get_diag_ind, i_vector_type&, 0, (empty_arg__));                    
    WRAPPER_METHOD_R_CONST (get_diag_ind_const, const i_vector_type&, 0, (empty_arg__));        
    WRAPPER_METHOD_R       (get_rows_ptr, i_vector_type&, 0, (empty_arg__));                    
    WRAPPER_METHOD_R_CONST (get_rows_ptr_const, const i_vector_type&, 0, (empty_arg__));        
    WRAPPER_METHOD_R       (get_cols_ind, i_vector_type&, 0, (empty_arg__));                    
    WRAPPER_METHOD_R_CONST (get_cols_ind_const, const i_vector_type&, 0, (empty_arg__));        
    WRAPPER_METHOD_R       (get_values, fp_storage_vector_type&, 0, (empty_arg__));             
    WRAPPER_METHOD_R_CONST (get_values_const, const fp_storage_vector_type&, 0, (empty_arg__)); 
    WRAPPER_METHOD_R_CONST (internal_check, int, 0, (empty_arg__));                             
  };
#endif //0
  PY_EXPORTER (py_bcsr_matrix_iface_exporter, py_matrix_iface_exporter)
    .def ("init_by_matrix",                     
          &T::init_by_matrix,
          args ("matrix"), "Initialize matrix by matrix :)")
    .def ("init",                               
          &T::init,
          args ("n_rows", "n_cols", "n_block_size", "n_non_zeros"), "Initialize matrix")
    .def ("init_struct",                        
          &T::init_struct,
          args ("n_rows", "n_cols", "n_non_zeros"), "Initialize matrix structure with out values")
    .def ("alloc_rows_ptr",                     
          &T::alloc_rows_ptr,
          args ("n_rows"), "Initialize rows_ptr vector")
    .def ("alloc_cols_ind",                     
          &T::alloc_cols_ind,
          args ("n_non_zeros"), "Initialize cols_ind vector")
    .def ("alloc_values",                       
          &T::alloc_values,
          args ("n_non_zeros", "n_block_size"), "Initialize values vector")
    .def ("alloc_cols_ind_and_values",          
          &T::alloc_cols_ind_and_values,
          args ("n_non_zeros", "n_block_size"), "Initialize values and cols_ind vectors")
    .def ("copy",                               
          &T::copy,
          args ("matrix"), "Initialize matrix by matrix and copy all content")
    .def ("get_rows_ptr",                       
          &T::get_rows_ptr, 
          args (""), "Return reference to the rows_ptr vector")
    .def ("get_cols_ind",                       
          &T::get_cols_ind, 
          args (""), "Return reference to the cols_ind vector")
    .def ("get_values",                         
          &T::get_values, 
          args (""), "Return reference to values vector")
    .def ("internal_check",                     
          &T::internal_check,
          args (""), "Return 0 if OK")
  PY_EXPORTER_END;

  PY_EXPORTER (py_bcsr_amg_matrix_iface_exporter, py_bcsr_matrix_iface_exporter)
    .def ("build_transpose",                    
          &T::build_transpose,
          args ("matrix", "rows_offset", "cols_offset", "new_n_rows"), "Build transpose matrix")
    .def ("build_transpose_struct",             
          &T::build_transpose_struct,
          args ("matrix", "rows_offset", "cols_offset", "new_n_rows"), "Build transpose matrix (only structure)")
    .def ("triple_matrix_product",             
          &T::triple_matrix_product,
          args ("R", "A", "P", "update_flag"), "calculate M = RAP (if update_flag == true calculate only values do not rebuild structure)")
  PY_EXPORTER_END;                               

  PY_EXPORTER (py_matrix_bcsr_tools_exporter, default_exporter)
    .def ("ascii_read_from_csr_format",         
          &T::ascii_read_from_csr_format,
          args ("matrix", "file"), "Read given matrix from ascii format")
    .def ("random_init",                        
          &T::random_init,
          args ("matrix", "n_rows", "n_block_size", "value_dispertion", "elems_in_row"),
          "Initialize matrix by random values")
    .def ("dense_init",                         
          &T::dense_init,
          args ("matrix", "n_rows", "n_block_size", "value_dispertion"),
          "Initialize matrix by random values (make dense matrix")
                                                 
  PY_EXPORTER_END;                               


  //PY_EXPORTER (py_dummy_exporter, default_exporter)
                                                 
  //PY_EXPORTER_END;                               
                                                 
  //! export matrices to python                  
  void py_export_bcsr_matrices ();                    

    
  } // namespace python
} // namespace blue_sky
#endif //
#endif //__PY_BCSR_H

