/**
 * \file matrix_vector_product.h
 * \brief matrix vector product impl for bcsr_matrix
 * \author Sergey Miryanov
 * \date 30.03.2009
 * */
#ifndef BS_MATRIX_MATRIX_VECTOR_PRODUCT_H_
#define BS_MATRIX_MATRIX_VECTOR_PRODUCT_H_

namespace blue_sky {

  // mul_alpha and alpha - get param from calc_lin_comb and reduce multiple loop in calc_lin_comb
  template <class item_array_t, class index_array_t, class v_vector_t, class r_vector_t, bool mul_alpha>
  int
  bcsr_matrix_vector_product_internal (const bcsr_matrix <item_array_t, index_array_t> *mx, const v_vector_t &v, r_vector_t &r, const typename r_vector_t::value_type &alpha = 1.0)
  {
    typedef typename index_array_t::value_type  i_type;
    typedef typename item_array_t::value_type   fp_type;
    typedef i_type                              index_t;
    typedef fp_type                             item_t;
    typedef typename v_vector_t::value_type     v_item_t;
    typedef typename r_vector_t::value_type     r_item_t;

    i_type i, j;
    i_type j1, j2;
    i_type cl               = 0;
    r_item_t *r_block       = 0;
    const fp_type *m_block  = 0;
    const v_item_t *v_block = 0;
    int b_sqr               = mx->n_block_size * mx->n_block_size;

    BS_ASSERT (v.size ());
    BS_ASSERT (r.size ());

    ////////////////////////////////////////////
    // loop through rows
    ////////////////////////////////////////////

    const index_t *rows = &mx->get_rows_ptr ()[0];
    const index_t *cols = &mx->get_cols_ind ()[0];
    const item_t *val   = &mx->get_values ()[0];
    const v_item_t *v_val = &v[0];
    r_item_t *r_val       = &r[0];

    if ( (mx->n_block_size) == 1)
      {
        r_item_t *r_block = &r_val[0];
        for (i = 0; i < (mx->n_rows); ++i, r_block += 1)
          {
            j1 = rows[i];
            j2 = rows[i + 1];

            ////////////////////////////////////////////
            // loop through columns in row
            ////////////////////////////////////////////
            const item_t *m_block = &val[j1];
            for (j = j1; j < j2; ++j, m_block += 1)
              {
                BS_ASSERT (j < (index_t)(mx->get_cols_ind()).size ());
                cl = cols[j];

                BS_ASSERT (cl * (mx->n_block_size) < (int)v.size ());
                v_block = &v_val[cl * 1];
                MV_PROD_1x1 (m_block, v_block, r_block);
              }
            if (mul_alpha)
              {
                r_block[0] *= alpha;
              }
          }
      }
    else if ( (mx->n_block_size) == 2)
      {
        for (i = 0; i < (mx->n_rows); ++i)
          {
            BS_ASSERT (i * (mx->n_block_size) < (index_t)r.size ());

            r_block = &r_val[i * (mx->n_block_size) ];
            j1      = rows[i];
            j2      = rows[i + 1];

            ////////////////////////////////////////////
            // loop through columns in row
            ////////////////////////////////////////////
            for (j = j1; j < j2; ++j)
              {
                BS_ASSERT (j < (index_t)(mx->get_cols_ind()).size ());
                cl = cols[j];
                m_block = &val[j * b_sqr];

                BS_ASSERT (cl * (mx->n_block_size) < (index_t)v.size ());
                v_block = &v_val[cl * (mx->n_block_size)];
                MV_PROD_2x2 (m_block, v_block, r_block);
              }
            if (mul_alpha)
              {
                r_block[0] *= alpha;
                r_block[1] *= alpha;
              }
          }
      }
    else if ( (mx->n_block_size) == 3)
      {
        r_item_t *r_block = &r_val[0];
        for (i = 0; i < (mx->n_rows); ++i, r_block += 3)
          {
            j1 = rows[i];
            j2 = rows[i + 1];

            ////////////////////////////////////////////
            // loop through columns in row
            ////////////////////////////////////////////
            const item_t *m_block = &val[j1 * 9];
            for (j = j1; j < j2; ++j, m_block += 9)
              {
                BS_ASSERT (j < (index_t)(mx->get_cols_ind()).size ());
                cl = cols[j];

                BS_ASSERT (cl * (mx->n_block_size) < (index_t)v.size ());
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

#ifdef CSR_MATRIX_VECTOR_PRODUCT_PARALLEL
#pragma omp parallel for private (r_block, j1, j2, j, cl, m_block, v_block)
#endif //CSR_MATRIX_VECTOR_PRODUCT_PARALLEL
    return 0;
  }

  template <class item_array_t, class index_array_t, class v_vector_t, class r_vector_t>
  int
  bcsr_matrix_vector_product (const bcsr_matrix <item_array_t, index_array_t> *mx, const v_vector_t &v, r_vector_t &r)
  {
    return bcsr_matrix_vector_product_internal <item_array_t, index_array_t, v_vector_t, r_vector_t, false> (mx, v, r);
  }

  template <class item_array_t, class index_array_t, class v_vector_t, class r_vector_t>
  int
  bcsr_matrix_vector_product_t (const bcsr_matrix <item_array_t, index_array_t> *mx, const v_vector_t &v, r_vector_t &r)
  {
    typedef typename index_array_t::value_type  i_type;
    typedef typename index_array_t::value_type  index_t;
    typedef typename item_array_t::value_type   item_t;

    i_type i, j1, j2, j, cl;
    const index_t *rows_ptr = &mx->get_rows_ptr ()[0];
    const index_t *cols_ind = &mx->get_cols_ind ()[0];
    const item_t *values    = &mx->get_values ()[0];

    BS_ASSERT (v.size ());
    BS_ASSERT (r.size ());
    BS_ASSERT (v.size () >= r.size ()) (v.size ()) (r.size ());
    BS_ASSERT (mx->n_rows == (int)v.size ());

    BS_ASSERT (mx->n_block_size <= 1);
    if (mx->n_block_size > 1)
      return -1;

    ////////////////////////////////////////////
    // loop through columns in transform matrix
    ////////////////////////////////////////////
    for (i = 0; i < mx->n_rows; ++i)
      {
        j1 = rows_ptr[i];
        j2 = rows_ptr[i + 1];

        ////////////////////////////////////////////
        // loop through rows in transform matrix
        ////////////////////////////////////////////
        for (j = j1; j < j2; ++j)
          {
            cl = cols_ind[j];
            r[cl] += v[i] * values[j];
          }
      }
    return 0;
  }

} // namespace blue_sky


#endif // #ifndef BS_MATRIX_MATRIX_VECTOR_PRODUCT_H_

