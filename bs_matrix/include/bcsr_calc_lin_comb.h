/**
 * \file calc_lin_comb.h
 * \brief calculate linear combination impl for bcsr matrix
 * \author Sergey Miryanov
 * \date 30.03.2009
 * */
#ifndef BS_MATRIX_CALC_LIN_COMB_H_
#define BS_MATRIX_CALC_LIN_COMB_H_

namespace blue_sky {

  template <typename item_array_t, typename index_array_t, typename rhs_array_t, typename sol_array_t>
  static int
  bcsr_calc_lin_comb (const bcsr_matrix <item_array_t, index_array_t> *matrix,
                 typename item_array_t::value_type          alpha,
                 typename item_array_t::value_type          beta,
                 const sol_array_t                          &u,
                 const rhs_array_t                          &v,
                 sol_array_t                                &r)
  {
    typedef typename sol_array_t::value_type    item_t;
    typedef typename index_array_t::value_type  index_t;

    static const item_t eps = item_t (1.0e-12);
    index_t i, n;
    int r_code = 0;

    item_t d = beta;
    if (fabs (alpha) > eps)
      {
        BS_ASSERT (u.size ());
        //BS_ASSERT (u.size () == r.size ());
        BS_ASSERT (r.size () >= (size_t)(matrix->n_rows * matrix->n_block_size));
        d /= alpha;
      }
    if (fabs (beta) > eps)
      {
        BS_ASSERT (v.size ());
        BS_ASSERT (v.size () == r.size ())(v.size ())(r.size ());
        BS_ASSERT (r.size () >= (size_t)(matrix->n_rows * matrix->n_block_size)) (r.size ()) (matrix->n_rows) (matrix->n_block_size) (matrix->n_rows * matrix->n_block_size);
      }

    if (fabs (beta) > eps)
      {
        i = 0;
        n = (index_t)r.size ();
        index_t n2 = n - (n % 4);
        // TODO:
#ifdef CSR_CALC_LIN_COMB_PARALLEL
#pragma omp parallel for
#endif
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
        assign (r, 0);
      }

    if (fabs (alpha) > eps)
      {
        //r_code = matrix->matrix_vector_product (u, r);
        bcsr_matrix_vector_product_internal <item_array_t, index_array_t, sol_array_t, sol_array_t, true> (matrix, u, r, alpha);
//        if (alpha != 1.0)
//          {
//            i = 0;
//            n = (index_t)r.size ();
//            index_t n2 = n - (n % 4);
//
//            // TODO:
//#ifdef CSR_CALC_LIN_COMB_PARALLEL
//#pragma omp parallel for
//#endif
//            for (; i < n2; i+=4)
//              {
//                r[i + 0] *= alpha;
//                r[i + 1] *= alpha;
//                r[i + 2] *= alpha;
//                r[i + 3] *= alpha;
//              }
//
//            for (; i < n; ++i)
//              r[i] *= alpha;
//          }
      }

    return r_code;
  }

} // namespace blue_sky

#endif // #ifndef BS_MATRIX_CALC_LIN_COMB_H_

