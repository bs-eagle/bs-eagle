#ifndef __MATRIX_IFACE_H
#define __MATRIX_IFACE_H
/** 
 * @file matrix_iface.h
 * @brief BS base interface class for matrix storage and manipulation
 * @author Oleg Borschuk
 * @date 2009-07-21
 */
#include <string>

#include "bs_assert.h"
#include "bs_tree.h"
#include "smart_ptr.h"
//#include "bsvector.h"
#include "bs_array.h"


namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  template <class strat_t>
  class matrix_iface: public bs_node
    {
    public:
      typedef typename strat_t::fp_type_t               fp_type_t;          //!< internal floating point type
      typedef typename strat_t::i_type_t                i_type_t;           //!< internal integer type
      //typedef typename strat_t::fp_vector_type          fp_vector_type_t;
      //typedef typename strat_t::i_vector_type           i_vector_type_t;

      //typedef ndarray<fp_type_t>                                fp_numpy_t;
      //typedef ndarray<fp_storage_type_t>                        fp_storage_numpy_t;
      //typedef ndarray<i_type_t>                                 i_numpy_t;
      typedef smart_ptr<bs_array<fp_type_t>, true>              sp_fp_array_t;
      typedef smart_ptr<bs_array<i_type_t>, true>               sp_i_array_t;

      //typedef smart_ptr<fp_storage_numpy_t, true>               sp_fp_storage_numpy_t;
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      
      //blue-sky class declaration
      //BLUE_SKY_TYPE_DECL_T(matrix_iface);


    public:
      //! destructor
      virtual ~matrix_iface ()
        {};

      //! calculate matrix vector product, v -- input vector, r -- output vector
      //! r += A * v
      //! return 0 if success
      virtual int matrix_vector_product (sp_fp_array_t v, sp_fp_array_t r) const = 0;
      
      //! calculate matrix^t vector product, v -- input vector, r -- output vector
      //! r += A^T * v
      //! return 0 if success
      virtual int matrix_vector_product_t (sp_fp_array_t v, sp_fp_array_t r) const = 0;

      //! calculate linear combination r = alpha * Au + beta * v
      //! alpha, beta -- scalar
      //! v, u -- input vector
      //! r -- output vector
      //! return 0 if success
      virtual int calc_lin_comb (fp_type_t alpha, fp_type_t beta, sp_fp_array_t u, sp_fp_array_t v, sp_fp_array_t r) const = 0;

      //! return total amount of allocated memory in bytes
      virtual fp_type_t get_allocated_memory_in_mbytes () const = 0;

      //! return block size 
      virtual i_type_t get_n_block_size () const = 0;

      //! return number of rows in matrix 
      virtual i_type_t get_n_rows () const = 0;

      //! return number of cols in matrix (return -1 in number of columns is unknown)
      virtual i_type_t get_n_cols () const = 0;

      //! return true if number of rows is equal to the number of columns
      virtual bool is_square () const = 0;
      
      //! initialize vector
      virtual void init_vector (sp_fp_array_t v) const = 0;
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN
    };

}//namespace blue_sky
#endif //__MATRIX_IFACE_H

