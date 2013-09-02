#ifndef __DENS_MATRIX_H
#define __DENS_MATRIX_H

#include "dens_matrix_iface.h"

namespace blue_sky
{
  /**
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN dens_matrix: public dens_matrix_iface
    {

    public:
      typedef matrix_iface                                      base_t;
      typedef smart_ptr<base_t, true>                           sp_base_t;
      typedef dens_matrix_iface                                 dens_matrix_iface_t;
      typedef smart_ptr<dens_matrix_iface_t, true>              sp_dens_matrix_iface_t;
      typedef dens_matrix                                       this_t;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (dens_matrix);
    public:

      //! destructor
      virtual ~dens_matrix ()
        {};

      //-----------------------------------------
      //  matrix_iface METHODS
      //-----------------------------------------

      //! calculate matrix vector product, v -- input vector, r -- output vector
      //! r += A * v
      //! return 0 if success
      virtual int matrix_vector_product (spv_double v, spv_double r) const;

      //! calculate matrix^t vector product, v -- input vector, r -- output vector
      //! r += A^T * v
      //! return 0 if success
      virtual int matrix_vector_product_t (spv_double v, spv_double r) const;

      //! calculate linear combination r = alpha * Au + beta * v
      //! alpha, beta -- scalar
      //! v, u -- input vector
      //! r -- output vector
      //! return 0 if success
      virtual int calc_lin_comb (t_double alpha, t_double beta, spv_double u, spv_double v, spv_double r) const;

      //! return total amount of allocated memory in bytes
      virtual t_double get_allocated_memory_in_mbytes () const
        {
          t_double d = 0;

          d += sizeof (this);
          d += sizeof (t_float) * values->size ();
          d /= 1024 * 1024;
          return d;
        }

      //! return block size
      virtual t_long get_n_block_size () const
      {
        return 1;
      }

      //! return number of rows in matrix
      virtual t_long get_n_rows () const
      {
        return n_rows;
      }

      //! return number of cols in matrix (return -1 in number of columns is unknown)
      virtual t_long get_n_cols () const
      {
        return n_cols;
      }

      //! return true if number of rows is equal to the number of columns
      virtual bool is_square () const
        {
          return n_cols == n_rows;
        }
      //! initialize vector
      virtual void init_vector (spv_double v) const
        {
          v->resize (n_rows);
          memset (&(*v)[0], 0, sizeof (t_double) * n_rows);
        }

      //-----------------------------------------
      //  dens_matrix_iface METHODS
      //-----------------------------------------

      /**
       * @brief initialize matrix
       *
       * @param matrix  -- <INPUT> referense to the matrix class
       *
       * @return 0 if success
       */
      virtual int init_by_matrix (sp_base_t matrix);

      /**
       * @brief initialize matrix
       *
       * @param new_n_rows              -- <INPUT> number of rows in matrix
       * @param new_n_cols              -- <INPUT> number of columns in matrix
       * @param calc_block_size         -- <INPUT> size of block used in matrix vector product (60 is a good choice)
       *
       * @return 0 if success
       */
      virtual int init (const t_long new_n_rows, const t_long new_n_cols, const t_long calc_block_size);

      /**
       * @brief make a copy of given matrix
       *
       * @param matrix  -- <INPUT> given matrix
       *
       * @return 0 if success
       */
      virtual int copy (sp_dens_matrix_iface_t matrix);

      /**
       * @brief return reference to the vector of values
       */
      virtual spv_float get_values ()
        {
          return values;
        }

      /**
       * @brief return const reference to the vector of values
       */
      //virtual const fp_storage_vector_type_t &get_values_const () const
      //  {
      //    return values;
      //  }


      /**
       * @brief return size of the calculation block in matrix vector product
       *        -1 -- use sequential algorithm
       */
      virtual t_long get_calc_block_size () const
        {
          return calc_block_size;
        }

      /**
       * @brief set size of the calculation block in matrix vector product
       *        if block_size * block_size can be put into the cache 1 memory
       *        this speedup calculation.
       *        60-80 is a good choice.
       *
       * @param block_size      -- <INPUT> new calculation block size (-1 for sequential calculation)
       */
      virtual void set_calc_block_size (const t_long block_size)
        {
          calc_block_size = block_size;
        }
      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      //! check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const
        {
          return 0;
        }

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const;
#endif //BSPY_EXPORTING_PLUGIN
    protected:
      void block_mv_product (t_long row_block, t_long col_block, const t_double *v, t_double *r) const;
      void block_mv_product_t (t_long row_block, t_long col_block, const t_double *v, t_double *r) const;

    protected:
      t_long n_rows;                  //!< number of rows in matrix
      t_long n_cols;                  //!< number of columns in matrix
      t_long calc_block_size;         //!< block size used in calculation method (don't affect on matrix storaging), -1 for sequential calculation
      //fp_storage_vector_type_t values;    //!< matrix values stored as linear vector by rows (i, j) -> (i * n_rows + j)
      spv_float values;
    };

}//namespace blue_sky
#endif //__DENS_MATRIX_H
