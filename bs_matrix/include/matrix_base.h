/**
 * @file matrix_base.h
 * @brief
 * @author Borschuk Oleg
 * @date 2008-03-25
 */
#ifndef MATRIX_BASE_H_
#define MATRIX_BASE_H_

#include "bs_object_base.h"
#include "bs_assert.h"
#include "throw_exception.h"


/**
 * @brief base class for matrix storage
 * <float, int>
 * <float, long>
 * <double, int>
 * <double, long>
 */


namespace blue_sky
{
  template <class fp_vector_type, class i_vector_type>
  class BS_API_PLUGIN matrix_base: public objbase
    {

    public:
      typedef typename fp_vector_type::value_type   fp_type_t;
      typedef typename i_vector_type::value_type    i_type_t;
      typedef fp_vector_type                        vec_fp_type;
      typedef fp_vector_type                        fp_vector_t;
      typedef i_vector_type                         i_vector_t;

      typedef fp_vector_type                        item_array_t;
      typedef i_vector_type                         index_array_t;
      typedef typename item_array_t::value_type     item_t;
      typedef typename index_array_t::value_type    index_t;

      typedef typename item_array_t::template array <float>::type   float_array_t;
      typedef typename item_array_t::template array <double>::type  double_array_t;

      typedef typename item_array_t::iterator       vec_iterator;
      typedef typename item_array_t::const_iterator vec_citerator;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL_T(matrix_base);


    public:
      //! destructor
      virtual ~matrix_base ()
      {
        n_block_size = 1;
        n_rows = 0;
        n_cols = 0;
      }

      // calculate matrix vector product, v -- input vector, r -- output vector
      // r += A * v
      // return 0 if success
      virtual int matrix_vector_product (const double_array_t & /*v*/, float_array_t & /*r*/) const
      {
        throw bs_exception ("matrix_vector_product", "PURE CALL");      
      }
      virtual int matrix_vector_product (const float_array_t & /*v*/, float_array_t & /*r*/) const
      {
        throw bs_exception ("matrix_vector_product", "PURE CALL");      
      }
      virtual int matrix_vector_product (const double_array_t & /*v*/, double_array_t & /*r*/) const
      {
        throw bs_exception ("matrix_vector_product", "PURE CALL");      
      }

      // calculate matrix^t vector product, v -- input vector, r -- output vector
      // r += A^t * v
      // return 0 if success
      virtual int matrix_vector_product_t (const double_array_t & /*v*/, float_array_t & /*r*/)
      {
        throw bs_exception ("matrix_vector_product_t", "PURE CALL");      
      }
      virtual int matrix_vector_product_t (const float_array_t & /*v*/, float_array_t & /*r*/)
      {
        throw bs_exception ("matrix_vector_product_t", "PURE CALL");      
      }
      virtual int matrix_vector_product_t (const double_array_t & /*v*/, double_array_t & /*r*/)
      {
        throw bs_exception ("matrix_vector_product_t", "PURE CALL");      
      }

      // calculate linear combination r = alpha * Au + beta * v
      virtual int calc_lin_comb (fp_type_t /*alpha*/, fp_type_t /*beta*/, const double_array_t & /*u*/, const float_array_t & /*v*/, double_array_t & /*r*/) const
      {
        throw bs_exception ("calc_lin_comb", "PURE CALL");      
      }
      virtual int calc_lin_comb (fp_type_t /*alpha*/, fp_type_t /*beta*/, const float_array_t & /*u*/, const double_array_t & /*v*/, float_array_t & /*r*/) const
      {
        throw bs_exception ("calc_lin_comb", "PURE CALL");      
      }
      virtual int calc_lin_comb (fp_type_t /*alpha*/, fp_type_t /*beta*/, const float_array_t & /*u*/, const float_array_t & /*v*/, float_array_t & /*r*/) const
      {
        throw bs_exception ("calc_lin_comb", "PURE CALL");      
      }
      virtual int calc_lin_comb (fp_type_t /*alpha*/, fp_type_t /*beta*/, const double_array_t & /*u*/, const double_array_t & /*v*/, double_array_t & /*r*/) const
      {
        throw bs_exception ("calc_lin_comb", "PURE CALL");      
      }

      // return total amount of allocated memory
      virtual i_type_t get_allocated_memory_in_bytes ()
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

      virtual index_array_t &
      get_rows_ptr ()
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual const index_array_t &
      get_rows_ptr () const
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual index_array_t &
      get_cols_ind ()
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual const index_array_t &
      get_cols_ind () const
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual item_array_t &
      get_values ()
      {
        static item_array_t dummy;
        return dummy;
      }

      virtual const item_array_t &
      get_values () const
      {
        static item_array_t dummy;
        return dummy;
      }

      i_type_t
      get_n_block_size () const
      {
        return n_block_size;
      }
      i_type_t
      get_n_rows () const 
      {
        return n_rows;
      }
      i_type_t
      get_n_cols () const
      {
        return n_cols;
      }

      void 
      set_n_block_size (i_type_t nb) 
      {
        n_block_size = nb;
      }
      void 
      set_n_rows (i_type_t n)
      {
        n_rows = n;
      }
      void 
      set_n_cols (i_type_t n)
      {
        n_cols = n;
      }

      virtual void
      init_vector (float_array_t &v)
      {
        bs_throw_exception ("PURE CALL");
      }
      virtual void
      init_vector (double_array_t &v)
      {
        bs_throw_exception ("PURE CALL");
      }
      virtual void
      init_vector (i_vector_t &v)
      {
        bs_throw_exception ("PURE CALL");
      }

    public:
      i_type_t n_block_size;        //!< block size
      i_type_t n_rows;              //!< number of rows
      i_type_t n_cols;              //!< number of columns
    };

}//namespace blue_sky

#endif//MATRIX_BASE_H_
