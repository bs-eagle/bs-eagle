#ifndef PY_MATRIX_IFACE_H_
#define PY_MATRIX_IFACE_H_
#include "matrix_iface.h"
#include "bdiag_matrix.h"
#include "jac_matrix.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include "python/py_bs_object_base.h"
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
    
  /** 
   * @brief python wrapper for matrix_iface class
   */
  CLASS_WRAPPER_T (1, (strat_t), matrix_iface, py_matrix_iface)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <strat_t>                              matrix_t;               
                                                                                        
    typedef typename matrix_t::fp_type_t                        fp_type_t;              
    typedef typename matrix_t::i_type_t                         i_type_t;        
    typedef typename matrix_t::fp_vector_type_t                 fp_vector_type_t;
    typedef typename matrix_t::fp_storage_vector_type_t         fp_storage_vector_type_t;
    typedef typename matrix_t::i_vector_type_t                  i_vector_type_t;
                                                                                        
                                                                                        
    typedef matrix_iface <strat_t>                              wrapped_t;              
    typedef BOOST_PP_CAT (py_matrix_iface, _base) <strat_t>     base_t;                 
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (1, (strat_t), py_matrix_iface);
  
  public:
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAP_PURE_METHOD_R_CONST (calc_lin_comb, int, 5, (fp_type_t, fp_type_t, const fp_vector_type_t&, const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAP_PURE_METHOD_R_CONST (get_allocated_memory_in_mbytes, fp_type_t, 0, (empty_arg__));       
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_block_size, i_type_t, 0, (empty_arg__));            
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_rows, i_type_t, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (get_n_cols, i_type_t, 0, (empty_arg__));                  
                                                                                        
    WRAP_PURE_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_CONST   (init_vector, void, 1, (fp_vector_type_t &));                       
    WRAP_PURE_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                   
  };
 
    
  /** 
   * @brief python wrapper for bdiag_matrix_iface class
   */
  CLASS_WRAPPER_T (1, (strat_t), bdiag_matrix_iface, py_bdiag_matrix_iface)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <strat_t>                              matrix_t;               
    typedef bdiag_matrix_iface <strat_t>                          bdiag_matrix_t;         
                                                                                        
    typedef typename matrix_t::fp_type_t                        fp_type_t;              
    typedef typename matrix_t::fp_storage_type_t                fp_storage_type_t;              
    typedef typename matrix_t::i_type_t                         i_type_t;              
    typedef typename matrix_t::fp_vector_type_t                 fp_vector_type_t;              
    typedef typename matrix_t::i_vector_type_t                  i_vector_type_t;              
    typedef typename matrix_t::fp_storage_vector_type_t         fp_storage_vector_type_t;              
                                                                                        
                                                                                        
    typedef bdiag_matrix_iface <strat_t>         wrapped_t;              
    typedef BOOST_PP_CAT (py_bdiag_matrix_iface, _base) <strat_t>        base_t;                 
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (1, (strat_t), py_bdiag_matrix_iface);
  
  public:
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAP_PURE_METHOD_R_CONST (calc_lin_comb, int, 5, (fp_type_t, fp_type_t, const fp_vector_type_t&, const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAP_PURE_METHOD_R_CONST (get_allocated_memory_in_mbytes, fp_type_t, 0, (empty_arg__));       
    WRAP_PURE_METHOD_R_CONST (get_n_block_size, i_type_t, 0, (empty_arg__));            
    WRAP_PURE_METHOD_R_CONST (get_n_rows, i_type_t, 0, (empty_arg__));                  
    WRAP_PURE_METHOD_R_CONST (get_n_cols, i_type_t, 0, (empty_arg__));                  
    WRAP_PURE_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_CONST   (init_vector, void, 1, (fp_vector_type_t &));                       
    WRAP_PURE_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                   
                                                                                        
    WRAP_PURE_METHOD_R       (init_by_matrix, int, 1, (const bdiag_matrix_t&));         
    WRAP_PURE_METHOD_R       (init, int, 2, (const i_type_t, const i_type_t));          
    WRAP_PURE_METHOD_R       (copy, int, 1, (const bdiag_matrix_t&));                   
    WRAP_PURE_METHOD_R       (get_diag, fp_storage_vector_type_t&, 0, (empty_arg__));     
    WRAP_PURE_METHOD_R_CONST (get_diag_const, const fp_storage_vector_type_t&, 0, (empty_arg__));          
  };

  /** 
   * @brief python wrapper for bdiag_matrix class
   */
  CLASS_WRAPPER_T (1, (strat_t), bdiag_matrix, py_bdiag_matrix)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <strat_t>        matrix_t;               
    typedef bdiag_matrix_iface <strat_t>   bdiag_matrix_t;         
                                                                                        
    typedef typename matrix_t::fp_type_t                        fp_type_t;              
    typedef typename matrix_t::fp_storage_type_t                fp_storage_type_t;              
    typedef typename matrix_t::i_type_t                         i_type_t;              
    typedef typename matrix_t::fp_vector_type_t                 fp_vector_type_t;              
    typedef typename matrix_t::i_vector_type_t                  i_vector_type_t;              
    typedef typename matrix_t::fp_storage_vector_type_t         fp_storage_vector_type_t;              
                                                                                        
                                                                                        
    typedef bdiag_matrix <strat_t>         wrapped_t;              
    typedef BOOST_PP_CAT (py_bdiag_matrix, _base) <strat_t>        base_t;                 
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (1, (strat_t), py_bdiag_matrix);

  public:
    WRAPPER_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAPPER_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAPPER_METHOD_R_CONST (calc_lin_comb, int, 5, (fp_type_t, fp_type_t, const fp_vector_type_t&, const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAPPER_METHOD_R_CONST (get_allocated_memory_in_mbytes, fp_type_t, 0, (empty_arg__));       
    WRAPPER_METHOD_R_CONST (get_n_block_size, i_type_t, 0, (empty_arg__));            
    WRAPPER_METHOD_R_CONST (get_n_rows, i_type_t, 0, (empty_arg__));                  
    WRAPPER_METHOD_R_CONST (get_n_cols, i_type_t, 0, (empty_arg__));                  
    WRAPPER_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAPPER_METHOD_CONST   (init_vector, void, 1, (fp_vector_type_t &));                       
    WRAPPER_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                  
                                                                                     
    WRAPPER_METHOD_R       (init_by_matrix, int, 1, (const bdiag_matrix_t&));        
    WRAPPER_METHOD_R       (init, int, 2, (const i_type_t, const i_type_t));         
    WRAPPER_METHOD_R       (copy, int, 1, (const bdiag_matrix_t&));                  
    WRAPPER_METHOD_R       (get_diag, fp_storage_vector_type_t&, 0, (empty_arg__));    
    WRAPPER_METHOD_R_CONST (get_diag_const, const fp_storage_vector_type_t&, 0, (empty_arg__));          
  };

  /** 
   * @brief python wrapper for jac_matrix_iface class
   */
  CLASS_WRAPPER_T (1, (strat_t), jac_matrix_iface, py_jac_matrix_iface)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <strat_t>              matrix_t;               
    typedef bdiag_matrix_iface <strat_t>        bdiag_matrix_t;       
    typedef bcsr_matrix_iface <strat_t>         bcsr_matrix_t;         
    typedef jac_matrix_iface <strat_t>          jac_matrix_t;           
    typedef smart_ptr<bcsr_matrix_t, true>      sp_bcsr_matrix_iface_t;    
    typedef smart_ptr<bdiag_matrix_t, true>     sp_bdiag_matrix_iface_t;   
                                                                                        
    typedef typename matrix_t::fp_type_t                        fp_type_t;              
    typedef typename matrix_t::fp_storage_type_t                fp_storage_type_t;              
    typedef typename matrix_t::i_type_t                         i_type_t;              
    typedef typename matrix_t::fp_vector_type_t                 fp_vector_type_t;              
    typedef typename matrix_t::i_vector_type_t                  i_vector_type_t;              
    typedef typename matrix_t::fp_storage_vector_type_t         fp_storage_vector_type_t;              
                                                                                        
    typedef jac_matrix_iface <strat_t>         wrapped_t;                    
    typedef BOOST_PP_CAT (py_jac_matrix_iface, _base) <strat_t> base_t;      
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (1, (strat_t), py_jac_matrix_iface);

  public:
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type_t&, fp_vector_type_t&));         
    WRAP_PURE_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type_t&, fp_vector_type_t&));       
    WRAP_PURE_METHOD_R_CONST (calc_lin_comb, int, 5, (fp_type_t, fp_type_t, const fp_vector_type_t&, const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAP_PURE_METHOD_R_CONST (get_allocated_memory_in_mbytes, fp_type_t, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_R_CONST (get_n_block_size, i_type_t, 0, (empty_arg__));            
    WRAP_PURE_METHOD_R_CONST (get_n_rows, i_type_t, 0, (empty_arg__));                  
    WRAP_PURE_METHOD_R_CONST (get_n_cols, i_type_t, 0, (empty_arg__));                  
    WRAP_PURE_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_CONST   (init_vector, void, 1, (fp_vector_type_t &));                       
    WRAP_PURE_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                   
                                                                                        
    WRAP_PURE_METHOD_R       (get_flux_matrix, sp_bcsr_matrix_iface_t, 0, (empty_arg__));                       
    WRAP_PURE_METHOD_R       (get_facility_matrix, sp_bcsr_matrix_iface_t, 0, (empty_arg__));                   
    WRAP_PURE_METHOD_R       (get_accum_matrix, sp_bdiag_matrix_iface_t, 0, (empty_arg__));                     
  };

  /** 
   * @brief python wrapper for jac_matrix class
   */
  CLASS_WRAPPER_T (1, (strat_t), jac_matrix, py_jac_matrix)
  {                                                                                     
  public:                                                                               
    typedef matrix_iface <strat_t>              matrix_t;               
    typedef bdiag_matrix_iface <strat_t>    bdiag_matrix_t;       
    typedef bcsr_matrix_iface <strat_t>    bcsr_matrix_t;         
    typedef jac_matrix_iface <strat_t>    jac_matrix_t;           
    typedef smart_ptr<bcsr_matrix_t, true>                   sp_bcsr_matrix_iface_t;    
    typedef smart_ptr<bdiag_matrix_t, true>                  sp_bdiag_matrix_iface_t;   
                                                                                        
    typedef typename matrix_t::fp_type_t                        fp_type_t;              
    typedef typename matrix_t::fp_storage_type_t                fp_storage_type_t;              
    typedef typename matrix_t::i_type_t                         i_type_t;              
    typedef typename matrix_t::fp_vector_type_t                 fp_vector_type_t;              
    typedef typename matrix_t::i_vector_type_t                  i_vector_type_t;              
    typedef typename matrix_t::fp_storage_vector_type_t         fp_storage_vector_type_t;              
                                                                                        
    typedef jac_matrix <strat_t>                                wrapped_t;                    
    typedef BOOST_PP_CAT (py_jac_matrix, _base) <strat_t> base_t;      
                                                                                        
  public:                                                                               
    CLASS_WRAPPER_DECL_T (1, (strat_t), py_jac_matrix);

  public:
    WRAPPER_METHOD_R_CONST (matrix_vector_product, int, 2, (const fp_vector_type_t&, fp_vector_type_t&));         
    WRAPPER_METHOD_R_CONST (matrix_vector_product_t, int, 2, (const fp_vector_type_t&, fp_vector_type_t&));       
    WRAPPER_METHOD_R_CONST (calc_lin_comb, int, 5, (fp_type_t, fp_type_t, const fp_vector_type_t&, const fp_vector_type_t&, fp_vector_type_t&)); 
    WRAPPER_METHOD_R_CONST (get_allocated_memory_in_mbytes, fp_type_t, 0, (empty_arg__));                       
    WRAPPER_METHOD_R_CONST (get_n_block_size, i_type_t, 0, (empty_arg__));            
    WRAPPER_METHOD_R_CONST (get_n_rows, i_type_t, 0, (empty_arg__));                  
    WRAPPER_METHOD_R_CONST (get_n_cols, i_type_t, 0, (empty_arg__));                  
    WRAPPER_METHOD_R_CONST (is_square, bool, 0, (empty_arg__));                       
    WRAPPER_METHOD_CONST   (init_vector, void, 1, (fp_vector_type_t &));                       
    WRAP_PURE_METHOD_R_CONST (py_str, std::string, 0, (empty_arg__));                   
                                                                                        
    WRAPPER_METHOD_R       (get_flux_matrix, sp_bcsr_matrix_iface_t, 0, (empty_arg__));                       
    WRAPPER_METHOD_R       (get_facility_matrix, sp_bcsr_matrix_iface_t, 0, (empty_arg__));                   
    WRAPPER_METHOD_R       (get_accum_matrix, sp_bdiag_matrix_iface_t, 0, (empty_arg__));                     
  };

  PY_EXPORTER (py_matrix_iface_exporter, default_exporter)
    .add_property ("n_block_size",              
                   &T::get_n_block_size, 
                   "Size of block for storing block matricies")
    .add_property ("n_rows",                    
                   &T::get_n_rows, 
                   "Number of rows in matrix")
    .add_property ("n_cols",                    
                   &T::get_n_cols, 
                   "Number of columns in matrix")  
    .add_property ("allocated_memory_in_mbytes",
                   &T::get_allocated_memory_in_mbytes, 
                   "Total amount of allocated memory in MBytes")
    .def ("matrix_vector_product",              
          &T::matrix_vector_product, 
          args ("v", "r"), "Matrix vector product (r += Av)")
    .def ("matrix_vector_product_t",            
          &T::matrix_vector_product_t, 
          args ("v", "r"), "Transpose matrix vector product (r += A^(T)v)")
    .def ("calc_lin_comb",                      
          &T::calc_lin_comb, 
          args ("alpha", "beta", "u", "v", "r"), "r = alpha * Au + beta * v")
    .def ("is_square",                          
          &T::is_square, 
          args (""), "Return true if matrix is square")
    .def ("init_vector",                        
          &T::init_vector, 
          args ("v"), "Initialize vector (for non MPI: call v.resize (n_rows, 0.0)")
    .def ("__str__",                            
          &T::py_str)
                                                 
  PY_EXPORTER_END;                               
 
  PY_EXPORTER (py_bdiag_matrix_iface_exporter, py_matrix_iface_exporter)
    .def ("init_by_matrix", 
          &T::init_by_matrix, 
          args ("matrix"), "Initialize matrix by matrix :)")
    .def ("init",           
          &T::init, 
          args ("n_rows", "n_block_size"), "Initialize matrix")
    .def ("copy",           
          &T::copy, 
          args ("matrix"), "Initialize matrix by matrix and copy all content")
    .def ("get_diag",       
          &T::get_diag, 
          args (""), "Return reference to the diagonal vector")
                                                 
  PY_EXPORTER_END;                               

  PY_EXPORTER (py_jac_matrix_iface_exporter, py_matrix_iface_exporter)
    .def ("get_flux_matrix",                    
          &T::get_flux_matrix,
          args (""), "Return smart pointer to the FLUX matrix")
    .def ("get_accum_matrix",                   
          &T::get_accum_matrix,
          args (""), "Return smart pointer to the Accomulative matrix")
    .def ("get_facility_matrix",                
          &T::get_facility_matrix,
          args (""), "Return smart pointer to the Facility matrix")
                                                 
  PY_EXPORTER_END;                               

  PY_EXPORTER (py_dummy_exporter, default_exporter)
                                                 
  PY_EXPORTER_END;                               
                                                 
  //! export matrices to python                  
  void py_export_matrices ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_MATRIX_BASE_H_
