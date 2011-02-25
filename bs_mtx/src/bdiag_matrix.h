/** 
 * @file bdiag_matrix.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-14
 */
#ifndef _BDIAG_MATRIX_H
#define _BDIAG_MATRIX_H
#include <string>
#include <sstream>

#include "bdiag_matrix_iface.h"

#include "matrix_macroses.h"
//#include "bdiag_matrix_impl.h"

namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN bdiag_matrix: public bdiag_matrix_iface
    {

    public:
      typedef matrix_iface                                      base_t;
      typedef bdiag_matrix_iface                                bdiag_matrix_iface_t;
      typedef bdiag_matrix                                      this_t;
      typedef smart_ptr<bdiag_matrix_iface_t, true>             sp_bdiag_matrix_iface_t;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (bdiag_matrix);

      //! destructor
      virtual ~bdiag_matrix ()
      {
      }

      //-----------------------------------------
      //  matrix_iface METHODS
      //-----------------------------------------

      //! calculate matrix vector product, v -- input vector, r -- output vector
      //! r += A * v
      //! return 0 if success
      virtual int matrix_vector_product (spv_double v, spv_double r) const
        {
          return matrix_vector_product_internal <false> (v, r);
        }

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
        d += sizeof (t_float) * diag->size ();
        d /= 1024 * 1024;
        return d;
      }

      //! return block size 
      virtual t_long get_n_block_size () const
      {
        return n_block_size;
      }

      //! return number of rows in matrix 
      virtual t_long get_n_rows () const 
      {
        return n_rows;
      }

      //! return number of cols in matrix (return -1 in number of columns is unknown)
      virtual t_long get_n_cols () const
      {
        return n_rows;
      }

      //! return true if number of rows is equal to the number of columns
      virtual bool is_square () const
        {
          return true;
        }
      //! initialize vector
      virtual void init_vector (spv_double v) const
        {
          v->resize (get_n_rows () * get_n_block_size ());
          memset (&(*v)[0], 0, v->size () * sizeof (t_double));

        }

      //-----------------------------------------
      //  bcsr_matrix_iface METHODS
      //-----------------------------------------

      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr, values
      virtual int init (const t_long new_n_rows, const t_long new_n_block_size)
        {
          if (new_n_block_size < 1)
            return -1;

          n_block_size = new_n_block_size;
          n_rows = new_n_rows;

          diag->resize (n_rows * n_block_size * n_block_size);
          memset (&(*diag)[0], 0, diag->size () * sizeof (t_float));
          return 0;
        }

      //! initialize by matrix
      virtual int init_by_matrix (sp_bdiag_matrix_iface_t matrix)
        {
          n_rows = matrix->get_n_rows ();
          n_block_size = matrix->get_n_block_size ();
          diag = matrix->get_diag ()->clone ();
          return 0;
        }

      //! copy all data from given matrix
      virtual int copy (sp_bdiag_matrix_iface_t matrix)
        {
          return init_by_matrix (matrix);
        }

      // ------------------------------------------------
      // GET matrix private DATA methods
      // ------------------------------------------------

      virtual spv_float get_diag ()
        {
          return diag;
        }

      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      // check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const
        {
          return 0;
        }
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "Block Diagonal matrix:\n";
          s << "\tmatrix size (blocks)     : " << get_n_rows () << "\n";
          s << "\tblock size               : " << get_n_block_size () << "\n";
          s << "\tallocated memory (Mb)    : " << get_allocated_memory_in_mbytes () << "\n";
          s << "-------------------------------------------------------------\n";

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
    protected:

      template<bool mul_alpha>
      int matrix_vector_product_internal (spv_double v, spv_double r, const t_double &alpha = 1.0) const;
      
      //__________________________________________
      //  VARIABLES
      //  ________________________________________
    public:

    protected:
      spv_float diag;        //!< matrix main diagonal (n_block_size * n_block_size * n_rows)

      t_long                  n_block_size;           //!< size of matrix block
      t_long                  n_rows;                 //!< number of rows in matrix
    };

/** 
 * @brief if template parameter of method is set to FALSE do r += Av
 *        else r +=Av, r *= alpha
 * 
 * @param v             -- <INPUT> given vector
 * @param r             -- <INPUT/OUTPUT> result vector
 * @param alpha         -- <INPUT> multiplier
 * 
 * @return 0 if success
 */
template<bool mul_alpha> int
bdiag_matrix::matrix_vector_product_internal (spv_double v, 
                                              spv_double r, 
                                              const t_double &alpha) const
{
  t_long i;
  t_long b_sqr = n_block_size * n_block_size;

  //BS_ASSERT (v.size ());
  //BS_ASSERT (r.size ());

  ////////////////////////////////////////////
  // loop through rows
  ////////////////////////////////////////////

  t_float *val   = &(*diag)[0];
  t_double *v_val         = &(*v)[0];
  t_double *r_val         = &(*r)[0];

  if (n_block_size == 1)
    {
      for (i = 0; i < n_rows; ++i)
        {
          MV_PROD_1x1 (val + i, v_val + i, r_val + i);
          if (mul_alpha)
            {
              r_val[i] *= alpha;
            }
        }
    }
  else if (n_block_size == 2)
    {
      for (i = 0; i < n_rows; ++i)
        {
          MV_PROD_2x2 (val + i * b_sqr, v_val + i * n_block_size, r_val + i * n_block_size);
          if (mul_alpha)
            {
              r_val[i * n_block_size] *= alpha;
              r_val[i * n_block_size + 1] *= alpha;
            }
        }
    }
  else if (n_block_size == 3)
    {
      for (i = 0; i < n_rows; ++i)
        {
          MV_PROD_3x3 (val + i * b_sqr, v_val + i * n_block_size, r_val + i * n_block_size);
          if (mul_alpha)
            {
              r_val[i * n_block_size + 0] *= alpha;
              r_val[i * n_block_size + 1] *= alpha;
              r_val[i * n_block_size + 2] *= alpha;
            }
        }
    }
  return 0;
}
}//namespace blue_sky
#endif //_BDIAG_MATRIX_H

