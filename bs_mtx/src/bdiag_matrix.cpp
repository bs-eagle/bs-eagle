/** 
 * @file bdiag_matrix.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-18
 */

#include "bs_matrix_stdafx.h"
#include "bdiag_matrix.h"
#include "strategies.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  template <class strat_t>
  bdiag_matrix<strat_t>::bdiag_matrix (bs_type_ctor_param) 
        : bdiag_matrix_iface<strat_t> (),
        diag (BS_KERNEL.create_object (fp_storage_array_t::bs_type ()))
    {
    }
  template <class strat_t>
  bdiag_matrix <strat_t>::bdiag_matrix (const bdiag_matrix & /*src*/) :bs_refcounter (),
        diag (BS_KERNEL.create_object (fp_storage_array_t::bs_type ()))
     {
     }
  template <class strat_t> int
  bdiag_matrix<strat_t>::matrix_vector_product_t (sp_fp_array_t v_, sp_fp_array_t r_) const
    {
      i_type_t i;
      fp_type_t *v = &(*v_)[0];
      fp_type_t *r = &(*r_)[0];
      fp_storage_type_t *d = &(*diag)[0];
      
      //BS_ASSERT (v.size ());
      //BS_ASSERT (r.size ());
      //BS_ASSERT (v.size () >= r.size ()) (v.size ()) (r.size ());
      //BS_ASSERT (n_rows == (i_type_t)v.size ());

      //BS_ASSERT (n_block_size <= 1);
      if (n_block_size > 1)
        return -1;

      ////////////////////////////////////////////
      // loop through columns in transform matrix
      ////////////////////////////////////////////
      for (i = 0; i < n_rows; ++i)
        {
          r[i] += v[i] * d[i];
        }
      return 0;
    }

  template <class strat_t> int
  bdiag_matrix<strat_t>::calc_lin_comb (fp_type_t alpha, fp_type_t beta, sp_fp_array_t u_, sp_fp_array_t v_, sp_fp_array_t r_) const
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
      }

    if (fabs (alpha) > eps)
      {
        matrix_vector_product_internal <true> (u_, r_, alpha);
      }

    return r_code;
  }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(bdiag_matrix, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(bdiag_matrix, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bdiag_matrix<base_strategy_fif>), 1,  (bdiag_matrix_iface <base_strategy_fif> ), "bdiag_matrix_fif", "Block Diag Matrix class", "Realization of Block Diag Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bdiag_matrix<base_strategy_did>), 1,  (bdiag_matrix_iface <base_strategy_did> ), "bdiag_matrix_did", "Block Diag Matrix class", "Realization of Block Diag Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bdiag_matrix<base_strategy_dif>), 1,  (bdiag_matrix_iface <base_strategy_dif> ), "bdiag_matrix_dif", "Block Diag Matrix class", "Realization of Block Diag Matricies", false);

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bdiag_matrix<base_strategy_flf>), 1,  (bdiag_matrix_iface <base_strategy_flf> ), "bdiag_matrix_flf", "Block Diag Matrix class", "Realization of Block Diag Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bdiag_matrix<base_strategy_dld>), 1,  (bdiag_matrix_iface <base_strategy_dld> ), "bdiag_matrix_dld", "Block Diag Matrix class", "Realization of Block Diag Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bdiag_matrix<base_strategy_dlf>), 1,  (bdiag_matrix_iface <base_strategy_dlf> ), "bdiag_matrix_dlf", "Block Diag Matrix class", "Realization of Block Diag Matricies", false);
}  // blue_sky namespace

