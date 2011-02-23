/** 
 * @file jac_matrix.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-14
 */

#include "bs_mtx_stdafx.h"
#include "jac_matrix.h"

namespace blue_sky
  {

  //! constructor
  
  jac_matrix::jac_matrix (bs_type_ctor_param /*param*/)
      : sp_accum_matrix (BS_KERNEL.create_object (bdiag_matrix_t::bs_type ()))
      , sp_flux_matrix (BS_KERNEL.create_object (bcsr_matrix_t::bs_type ()))
      , sp_facility_matrix (BS_KERNEL.create_object (bcsr_matrix_t::bs_type ()))
  {

  }

  //! copy constructor
  
  jac_matrix::jac_matrix (const jac_matrix &matrix) 
        : bs_refcounter () 
  {
    BS_ASSERT (false && "TEST ME");

    if (this != &matrix)
      {
        *sp_accum_matrix = *(matrix.sp_accum_matrix);
        *sp_facility_matrix = *(matrix.sp_facility_matrix);
        *sp_flux_matrix = *(matrix.sp_flux_matrix);
      }

  }
  int
  jac_matrix::calc_lin_comb (t_double alpha, 
                             t_double beta, 
                             spv_double u_,
                             spv_double v_, 
                             spv_double r_) const
  {
    static const t_double eps = t_double (1.0e-12);
    t_long i, n;
    int r_code = 0;

    t_double *v = &(*v_)[0];
    t_double *r = &(*r_)[0];

    t_double d = beta;
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

    n = (t_long)r_->size ();
    if (fabs (beta) > eps)
      {
        i = 0;
        t_long n2 = n - (n % 4);
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
        memset (r, 0, sizeof (t_double) * n);
        //r.assign (n, 0);
      }

    if (fabs (alpha) > eps)
      {
        matrix_vector_product (u_, r_);
        i = 0;
        n = (t_long)r_->size ();
        t_long n2 = n - (n % 4);
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

  BLUE_SKY_TYPE_STD_CREATE (jac_matrix);
  BLUE_SKY_TYPE_STD_COPY (jac_matrix);


  BLUE_SKY_TYPE_IMPL (jac_matrix, jac_matrix_iface, "jac_matrix", "Jacobian Matrix class", "Realization of Jacobian Matricies");
  }


