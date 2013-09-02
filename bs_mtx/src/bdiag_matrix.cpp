/** 
 * @file bdiag_matrix.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-18
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "bdiag_matrix.h"
#include "strategies.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  bdiag_matrix::bdiag_matrix (bs_type_ctor_param) 
        : bdiag_matrix_iface (),
        diag (BS_KERNEL.create_object (v_float::bs_type ()))
    {
    }
  bdiag_matrix::bdiag_matrix (const bdiag_matrix & /*src*/) :bs_refcounter (),
        diag (BS_KERNEL.create_object (v_float::bs_type ()))
     {
     }
  int
  bdiag_matrix::matrix_vector_product_t (spv_double v_, spv_double r_) const
    {
      t_long i;
      t_double *v = &(*v_)[0];
      t_double *r = &(*r_)[0];
      t_float *d = &(*diag)[0];
      
      //BS_ASSERT (v.size ());
      //BS_ASSERT (r.size ());
      //BS_ASSERT (v.size () >= r.size ()) (v.size ()) (r.size ());
      //BS_ASSERT (n_rows == (t_long)v.size ());

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

  int
  bdiag_matrix::calc_lin_comb (t_double alpha, t_double beta, spv_double u_, spv_double v_, spv_double r_) const
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
      }

    if (fabs (alpha) > eps)
      {
        matrix_vector_product_internal <true> (u_, r_, alpha);
      }

    return r_code;
  }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (bdiag_matrix);
  BLUE_SKY_TYPE_STD_COPY (bdiag_matrix);

  BLUE_SKY_TYPE_IMPL (bdiag_matrix, bdiag_matrix_iface, "bdiag_matrix", "Block Diag Matrix class", "Realization of Block Diag Matricies");
}  // blue_sky namespace

