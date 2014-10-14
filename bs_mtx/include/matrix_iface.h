#ifndef __MATRIX_IFACE_H
#define __MATRIX_IFACE_H
/** 
 * @file matrix_iface.h
 * @brief BS base interface class for matrix storage and manipulation
 * @author Oleg Borschuk
 * @date 2009-07-21
 */
#include "bs_assert.h"
#include "bs_tree.h"
#include "smart_ptr.h"
#include "bs_array.h"
#include "conf.h"

#include <string>

namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  class BS_API_PLUGIN matrix_iface: public bs_node
    {
    public:

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      
    public:
      //! destructor
      virtual ~matrix_iface ()
        {};

      //! calculate matrix vector product, v -- input vector, r -- output vector
      //! r += A * v
      //! return 0 if success
      virtual int matrix_vector_product (spv_double v, spv_double r) const = 0;
      
      //! calculate matrix^t vector product, v -- input vector, r -- output vector
      //! r += A^T * v
      //! return 0 if success
      virtual int matrix_vector_product_t (spv_double v, spv_double r) const = 0;

      //! calculate linear combination r = alpha * Au + beta * v
      //! alpha, beta -- scalar
      //! v, u -- input vector
      //! r -- output vector
      //! return 0 if success
      virtual int calc_lin_comb (t_double alpha, t_double beta, spv_double u, spv_double v, spv_double r) const = 0;

      //! return total amount of allocated memory in bytes
      virtual t_double get_allocated_memory_in_mbytes () const = 0;

      //! return block size 
      virtual t_long get_n_block_size () const = 0;

      //! return number of rows in matrix 
      virtual t_long get_n_rows () const = 0;

      //! return number of cols in matrix (return -1 in number of columns is unknown)
      virtual t_long get_n_cols () const = 0;

      //! return true if number of rows is equal to the number of columns
      virtual bool is_square () const = 0;
      
      //! initialize vector
      virtual void init_vector (spv_double v) const = 0;
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN
    };

}//namespace blue_sky
#endif //__MATRIX_IFACE_H

