#ifndef _JAC_MATRIX_H
#define _JAC_MATRIX_H

#include <string>
#include <sstream>

#include "jac_matrix_iface.h"
#include "bdiag_matrix.h"
#include "bcsr.h"


namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  template <class strat_t>
  class BS_API_PLUGIN jac_matrix: public jac_matrix_iface<strat_t>
    {
    public:

      typedef typename strat_t::fp_vector_type                  fp_vector_type_t;
      typedef typename strat_t::i_vector_type                   i_vector_type_t;
      typedef typename strat_t::fp_storage_vector_type          fp_storage_vector_type_t;
      typedef typename strat_t::fp_type_t                       fp_type_t;
      typedef typename strat_t::i_type_t                        i_type_t;
      typedef typename strat_t::fp_storage_type_t               fp_storage_type_t;

      typedef bcsr_matrix_iface<strat_t>                        bcsr_matrix_iface_t;
      typedef bdiag_matrix_iface<strat_t>                       bdiag_matrix_iface_t;

      typedef bcsr<strat_t>                                     bcsr_matrix_t;
      typedef bdiag_matrix<strat_t>                             bdiag_matrix_t;

      typedef smart_ptr<bcsr_matrix_t, true>                    sp_bcsr_matrix_t;
      typedef smart_ptr<bdiag_matrix_t, true>                   sp_bdiag_matrix_t;

      typedef smart_ptr<bcsr_matrix_iface_t, true>              sp_bcsr_matrix_iface_t;
      typedef smart_ptr<bdiag_matrix_iface_t, true>             sp_bdiag_matrix_iface_t;

      typedef jac_matrix<strat_t>                               this_t;
      typedef matrix_iface<strat_t>                             base_t;

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;


      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL_T(jac_matrix);

      //! destructor
      virtual ~jac_matrix ()
      {
      }

      //-----------------------------------------
      //  matrix_iface METHODS
      //-----------------------------------------

      //! calculate matrix vector product, v -- input vector, r -- output vector
      //! r += A * v
      //! return 0 if success
      virtual int matrix_vector_product (sp_fp_array_t v, sp_fp_array_t r) const
      {
        return sp_accum_matrix->matrix_vector_product (v, r) 
               || sp_flux_matrix->matrix_vector_product (v, r) 
               || sp_facility_matrix->matrix_vector_product (v, r);
      }

      //! calculate matrix^t vector product, v -- input vector, r -- output vector
      //! r += A^T * v
      //! return 0 if success
      virtual int matrix_vector_product_t (sp_fp_array_t v, sp_fp_array_t r) const
      {
        return sp_accum_matrix->matrix_vector_product_t (v, r) 
               || sp_flux_matrix->matrix_vector_product_t (v, r) 
               || sp_facility_matrix->matrix_vector_product_t (v, r);
      }

      //! calculate linear combination r = alpha * Au + beta * v
      //! alpha, beta -- scalar
      //! v, u -- input vector
      //! r -- output vector
      //! return 0 if success
      virtual int calc_lin_comb (fp_type_t alpha, fp_type_t beta, sp_fp_array_t u, sp_fp_array_t v, sp_fp_array_t r) const;

      //! return total amount of allocated memory in bytes
      virtual fp_type_t get_allocated_memory_in_mbytes () const
      {
        return sp_accum_matrix->get_allocated_memory_in_mbytes ()
               + sp_flux_matrix->get_allocated_memory_in_mbytes ()
               + sp_facility_matrix->get_allocated_memory_in_mbytes ();
      }

      //! return block size 
      virtual i_type_t get_n_block_size () const
      {
        return sp_flux_matrix->get_n_block_size ();
      }

      //! return number of rows in matrix 
      virtual i_type_t get_n_rows () const 
      {
        return sp_flux_matrix->get_n_rows ();
      }

      //! return number of cols in matrix (return -1 in number of columns is unknown)
      virtual i_type_t get_n_cols () const
      {
        return sp_flux_matrix->get_n_cols ();
      }

      //! return true if number of rows is equal to the number of columns
      virtual bool is_square () const
        {
          return sp_flux_matrix->is_square ();
        }
      //! initialize vector
      virtual void init_vector (sp_fp_array_t v) const
        {
          v->resize (sp_flux_matrix->get_n_rows () * sp_flux_matrix->get_n_block_size ());
          memset (&(*v)[0], 0, sizeof (fp_type_t) * v->size ());
        }

      //-----------------------------------------
      //  jac_matrix_iface METHODS
      //-----------------------------------------

      // ------------------------------------------------
      // GET matrix private DATA methods
      // ------------------------------------------------

      //! return flux (regular) matrix
      virtual sp_bcsr_matrix_iface_t get_flux_matrix ()
        {
          return sp_flux_matrix;
        }

      //! return facility (irregular) matrix
      virtual sp_bcsr_matrix_iface_t get_facility_matrix ()
        {
          return sp_facility_matrix;
        }

      //! return accumulative matrix
      virtual sp_bdiag_matrix_iface_t get_accum_matrix ()
        {
          return sp_accum_matrix;
        }

      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      // check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const
        {
          return sp_accum_matrix->internal_check () 
                 || sp_flux_matrix->internal_check () 
                 || sp_facility_matrix->internal_check ();
        }
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "Accumulative part of jacobian matrix:\n";
          s << sp_accum_matrix->py_str ();
          s << "Flux part of jacobian matrix:\n";
          s << sp_flux_matrix->py_str ();
          s << "Facility part of jacobian matrix:\n";
          s << sp_facility_matrix->py_str ();

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
    public:

    protected:
      
      sp_bdiag_matrix_t sp_accum_matrix;
      sp_bcsr_matrix_t sp_flux_matrix;
      sp_bcsr_matrix_t sp_facility_matrix;

    };

template <class strat_t> int
jac_matrix<strat_t>::calc_lin_comb (fp_type_t alpha, 
                                    fp_type_t beta, 
                                    sp_fp_array_t u_,
                                    sp_fp_array_t v_, 
                                    sp_fp_array_t r_) const
{
  static const fp_type_t eps = fp_type_t (1.0e-12);
  i_type_t i, n;
  int r_code = 0;

  fp_type_t *v = &(*v_)[0];
  fp_type_t *r = &(*r_)[0];

  fp_type_t d = beta;
  if (fabs (alpha) > eps)
    {
      //BS_ASSERT (u.size ());
      //BS_ASSERT (r.size () >= (size_t)(n_rows * n_block_size));
      d /= alpha;
    }
  if (fabs (beta) > eps)
    {
      //BS_ASSERT (v.size ());
      //BS_ASSERT (v.size () == r.size ())(v.size ())(r.size ());
      //BS_ASSERT (r.size () >= (size_t)(n_rows * n_block_size)) (r.size ()) (n_rows) (n_block_size) (n_rows * n_block_size);
    }

  n = (i_type_t)r_->size ();
  if (fabs (beta) > eps)
    {
      i = 0;
      i_type_t n2 = n - (n % 4);
      // TODO:
      for (; i < n2; i+=4)
        {
          r[i + 0] = v[i + 0] * d;
          r[i + 1] = v[i + 1] * d;
          r[i + 2] = v[i + 2] * d;
          r[i + 3] = v[i + 3] * d;
        }

      for (; i < n; ++i)
        r[i] = v[i] * d;
    }
  else
    {
      memset (r, 0, sizeof (fp_type_t) * n);
      //r.assign (n, 0);
    }

  if (fabs (alpha) > eps)
    {
      matrix_vector_product (u_, r_);
      i = 0;
      n = (i_type_t)r_->size ();
      i_type_t n2 = n - (n % 4);
      // TODO:
      for (; i < n2; i+=4)
        {
          r[i + 0] *= alpha;
          r[i + 1] *= alpha;
          r[i + 2] *= alpha;
          r[i + 3] *= alpha;
        }

      for (; i < n; ++i)
        r[i] *= alpha;
    }

  return r_code;
}


}//namespace blue_sky

#endif //_JACOBIAN_H

