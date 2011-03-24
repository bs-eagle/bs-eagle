/** 
 * @file bcsr.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-13
 */
#ifndef __BCSR_H
#define __BCSR_H

#include "matrix_macroses.h"

#include "bcsr_amg_matrix_iface.h"
//#include "bcsr_matrix_impl.h"

#include <vector>
#include <string>

namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN bcsr: public bcsr_amg_matrix_iface
    {

    public:

      typedef matrix_iface                                      base_t;
      typedef bcsr_matrix_iface                                 bcsr_matrix_iface_t;
      typedef smart_ptr <bcsr_matrix_iface_t, true>             sp_bcsr_matrix_iface_t;
      typedef bcsr                                              this_t;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (bcsr);
    public:
      
      //! destructor
      virtual ~bcsr ()
        {};

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
      virtual int calc_lin_comb (t_double alpha, t_double beta, spv_double u, 
                                 spv_double v, spv_double r) const;

      //! return total amount of allocated memory in bytes
      virtual t_double get_allocated_memory_in_mbytes () const;

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
          //v->resize (get_n_rows () * get_n_block_size (), 0);
          v->resize (get_n_rows () * get_n_block_size ());
          memset (&(*v)[0], 0, sizeof (t_double) * v->size ());
        }

      //-----------------------------------------
      //  bcsr_matrix_iface METHODS
      //-----------------------------------------

      /** 
       * @brief initialize matrix by numpy arrays
       * 
       * @param n_rows_         -- <INPUT> number of rows
       * @param n_cols_         -- <INPUT> number of columns
       * @param n_block_size_   -- <INPUT> sizeof matrix block 
       * @param rows_           -- <INPUT> numpy array (rows_ptr information)
       * @param cols_           -- <INPUT> numpy array (cols_ind information)
       * @param values_         -- <INPUT> numpy array (values information)
       * @param make_copy       -- <INPUT> if true -- copy memory of rows_, cols_, values_; if false use the same memory
       * 
       * @return 0 if success
       */
      virtual int init_by_arrays (const t_long n_rows_, const t_long n_cols_, const t_long n_block_size_,
                                  spv_long rows_, spv_long cols_, spv_float values_, bool make_copy)
        {
          n_rows = n_rows_;
          n_cols = n_cols_;
          n_block_size = n_block_size_;
          if (make_copy)
            {
              rows_ptr = rows_->clone ();
              cols_ind = cols_->clone ();
              values = values->clone ();
            }
          else
            {
              rows_ptr = rows_;
              cols_ind = cols_;
              values = values_;
            }
          return 0;
        }

      /** 
       * @brief initialize matrix by raw pointers
       * 
       * @param n_rows_         -- <INPUT> number of rows
       * @param n_cols_         -- <INPUT> number of columns
       * @param n_block_size_   -- <INPUT> sizeof matrix block 
       * @param rows_           -- <INPUT> rows_ptr information
       * @param cols_           -- <INPUT> cols_ind information
       * @param values_         -- <INPUT> values information
       * @param make_copy       -- <INPUT> if true -- copy memory of rows_, cols_, values_; if false use the same memory
       * 
       * @return 0 if success
       */
      virtual int init_by_raw_ptr (const t_long n_rows_, const t_long n_cols_, const t_long n_block_size_,
                                  t_long *rows_, t_long *cols_, t_float *values_, bool make_copy)
        {
          if (init (n_rows_, n_cols_, n_block_size_, rows_[n_rows_]))
            return -1;
          
          if (make_copy)
            {
              int nnz = rows_[n_rows_];

              memcpy (&(*rows_ptr)[0], rows_, sizeof (t_long) * (n_rows_ + 1));
              memcpy (&(*cols_ind)[0], cols_, sizeof (t_long) * nnz);
              memcpy (&(*values)[0], values_, sizeof (t_float) * nnz * n_block_size_ * n_block_size_);
            }
          return 0;
        }

      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr, values
      virtual int init_by_matrix (sp_bcsr_matrix_iface_t m);

      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr, values
      virtual int init (const t_long new_n_rows, const t_long new_n_cols, const t_long new_n_blok_size,
                        const t_long new_n_non_zeros);

      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr
      virtual int init_struct (const t_long new_n_rows, const t_long new_n_cols, const t_long new_n_non_zeros);

      //! setup n_rows, and allocate memory for rows_ptr
      virtual int alloc_rows_ptr (const t_long new_n_rows);

      //! allocate memory for column indexes and set it with default values
      virtual int alloc_cols_ind (const t_long new_n_non_zeros);

      //! allocate memory for values and set it with default values
      virtual int alloc_values (const t_long new_n_non_zeros, const t_long new_n_blok_size);

      //! allocate memory for cols_ind and values
      virtual int alloc_cols_ind_and_values (const t_long new_n_non_zeros, const t_long new_n_blok_size);

      virtual void set_n_cols (const t_long new_n_cols)
        {
          n_cols = new_n_cols;
        }

      //! TODO should be removed
      virtual int copy (sp_bcsr_matrix_iface_t m)
        {
          return init_by_matrix (m);
        }

      //! build B = A^T from given csr matrix, return 0 if success,
      //! rows number and offset can be set manually with new_n_rows and rows_offset
      virtual int build_transpose (sp_bcsr_matrix_iface_t m, const t_long rows_offset/* = 0*/, 
                                   const t_long cols_offset/* = 0*/, const t_long new_n_rows/* = 0*/);

      //! build B = A^T from given csr matrix using matrix structure only, without values, return 0 if success,
      //! rows number and offset can be set manually with new_n_rows and rows_offset
      virtual int build_transpose_struct (sp_bcsr_matrix_iface_t m, const t_long rows_offset/* = 0*/,
                                          const t_long cols_offset/* = 0*/, const t_long new_n_rows/* = 0*/);
      /** 
       * @brief realize M = RAP (triple matrix product)
       * 
       * @param r_matrix        -- <INPUT> R matrix
       * @param a_matrix        -- <INPUT> A matrix
       * @param p_matrix        -- <INPUT> P matrix
       * @param update          -- <INPUT> if true than use existing matrix structure and update only values
       * 
       * @return 0 if success
       */
      virtual int triple_matrix_product (sp_bcsr_matrix_iface_t r_matrix, sp_bcsr_matrix_iface_t a_matrix, 
                                         sp_bcsr_matrix_iface_t p_matrix, const bool update);
      // ------------------------------------------------
      // GET matrix private DATA methods
      // ------------------------------------------------
      //! return number of nonzeros elements
      virtual t_long get_n_non_zeros () const
        {
          int l = rows_ptr->size (); 
          if (l)
            return (*rows_ptr)[l - 1];
          else
            return 0;
        }

      virtual spv_long get_rows_ptr ()
      {
        return rows_ptr;
      }

      //virtual const i_vector_type_t &get_rows_ptr_const () const
      //{
      //  return impl.get_rows_ptr_const ();
      //}

      virtual spv_long get_cols_ind ()
      {
        return cols_ind;
      }

      //virtual const i_vector_type_t &get_cols_ind_const () const
      //{
      //  return impl.get_cols_ind_const ();
      //}

      virtual spv_float get_values ()
      {
        return values;
      }

      //virtual const fp_storage_vector_type_t &get_values_const () const
      //{
      //  return impl.get_values_const ();
      //}

      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      // check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const
        {
          return 0; // TODO: implement
        }
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "Block Compressed Sparce Row matrix:\n";
          if (n_rows > 0)
            {
              s << "\tNumber of rows:             " << n_rows << "\n";
              s << "\tNumber of columns:          " << n_cols << "\n";
              s << "\tBlock size:                 " << n_block_size << "\n";
              s << "\tNumber of non zero blocks:  " << get_n_non_zeros () << "\n"; 
              s << "\tAllocated memory (MBytes):  " << get_allocated_memory_in_mbytes () << "\n";
            }
          else
            {
              s << "\tmatrix have not been initialized.\n";
            }
          s << "-----------------------------------------------------------\n";

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
    public:
    
    protected:
      template<bool mul_alpha>
      int matrix_vector_product_internal (spv_double v, spv_double r, const t_double &alpha = 1.0) const;

    //______________________________________
    //  VARIABLES
    //______________________________________
    protected:
      
      t_long                    n_block_size;           //!< size of matrix block
      t_long                    n_rows;                 //!< number of rows in matrix
      t_long                    n_cols;                 //!< number of columns in matrix
      spv_float                 values;                 //!< values array
      spv_long                  cols_ind;               //!< cols array
      spv_long                  rows_ptr;               //!< rows array

      // workspace
      std::vector<t_long> a_marker_vec;
      std::vector<t_long> p_marker_vec;

      //bcsr_matrix_impl<fp_vector_type_t, i_vector_type_t, fp_storage_vector_type_t> impl;    
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
template <bool mul_alpha> int
bcsr::matrix_vector_product_internal (spv_double v_, 
                                      spv_double r_, 
                                      const t_double &alpha) const
{
  t_long i, j;
  t_long j1, j2;
  t_long cl               = 0;
  t_double *r_block = 0;
  const t_float *m_block = 0;
  const t_double *v_block = 0;
  t_long b_sqr = n_block_size * n_block_size;

  t_long *rows = &(*rows_ptr)[0];
  t_long *cols = &(*cols_ind)[0];
  t_float *val = &(*values)[0];
  t_double *v_val = &(*v_)[0];
  t_double *r_val = &(*r_)[0];

  //BS_ASSERT (v.size ());
  //BS_ASSERT (r.size ());

  ////////////////////////////////////////////
  // loop through rows
  ////////////////////////////////////////////

  if (n_block_size == 1)
    {
      r_block = &r_val[0];
      for (i = 0; i < n_rows; ++i, r_block += 1)
        {
          j1 = rows[i];
          j2 = rows[i + 1];

          ////////////////////////////////////////////
          // loop through columns in row
          ////////////////////////////////////////////
          m_block = &val[j1];
          for (j = j1; j < j2; ++j, m_block += 1)
            {
              //BS_ASSERT (j < cols_ind.size ());
              cl = cols[j];

              //BS_ASSERT (cl * n_block_size < (t_long)v.size ());
              v_block = &v_val[cl * 1];
              MV_PROD_1x1 (m_block, v_block, r_block);
            }
          if (mul_alpha)
            {
              r_block[0] *= alpha;
            }
        }
    }
  else if (n_block_size == 2)
    {
      for (i = 0; i < n_rows; ++i)
        {
          //BS_ASSERT (i * n_block_size < (t_long)r.size ());

          r_block = &r_val[i * n_block_size];
          j1      = rows[i];
          j2      = rows[i + 1];

          ////////////////////////////////////////////
          // loop through columns in row
          ////////////////////////////////////////////
          for (j = j1; j < j2; ++j)
            {
              //BS_ASSERT (j < (t_long)cols_ind.size ());
              cl = cols[j];
              m_block = &val[j * b_sqr];

              //BS_ASSERT (cl * n_block_size < (t_long)v.size ());
              v_block = &v_val[cl * n_block_size];
              MV_PROD_2x2 (m_block, v_block, r_block);
            }
          if (mul_alpha)
            {
              r_block[0] *= alpha;
              r_block[1] *= alpha;
            }
        }
    }
  else if (n_block_size == 3)
    {
      r_block = &r_val[0];
      for (i = 0; i < n_rows; ++i, r_block += 3)
        {
          j1 = rows[i];
          j2 = rows[i + 1];

          ////////////////////////////////////////////
          // loop through columns in row
          ////////////////////////////////////////////
          m_block = &val[j1 * 9];
          for (j = j1; j < j2; ++j, m_block += 9)
            {
              //BS_ASSERT (j < (t_long)cols_ind.size ());
              cl = cols[j];

              //BS_ASSERT (cl * n_block_size < (t_long)v.size ());
              v_block = &v_val[cl * 3];
              MV_PROD_3x3 (m_block, v_block, r_block);
            }

          if (mul_alpha)
            {
              r_block[0] *= alpha;
              r_block[1] *= alpha;
              r_block[2] *= alpha;
            }
        }
    }
  return 0;
}
}//namespace blue_sky
#endif //__BCSR_H


