/**
 * @file mbcsr_matrix.cpp
 * @brief Block CSR multi matrix implementation
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-27
 */
#include "mbcsr_matrix.h"
#include "matrix_macroses.h"
#include <vector>
#include <memory.h>

namespace blue_sky
  {

  //! constructor

  mbcsr_matrix::mbcsr_matrix (bs_type_ctor_param /*param*/)
  {
  }

  //! copy constructor
  mbcsr_matrix::mbcsr_matrix (const mbcsr_matrix &matrix)
        : bs_refcounter ()
  {
    BS_ASSERT (false && "TEST ME");

    if (this != &matrix)
      {
        mat_map = matrix.mat_map;
      }
  }

  int
  mbcsr_matrix::calc_lin_comb (t_double alpha,
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

  /**
   * @brief merge matrices in mat_map (assume matrix is square and diagonal element exist)
   * @param filter -- input flags -- if (filter[i]) use row i
   * @return smart_ptr to merged matrix
   */

  mbcsr_matrix::sp_bcsr_iface_t
  mbcsr_matrix::merge (spv_int filter)
  {
    sp_bcsr_iface_t m;

    if (mat_map.size ())
      {
        map_t::iterator it, ie, ib;
        t_long nr, nb, b_sqr;
        std::vector<t_long> flg;
        t_long *rows, *i_rows;
        t_long *cols, *i_cols;
        t_float *values, *i_values;
        t_int *flt = 0;
        t_long j1, j2, j, cl;
        t_long count;

        if (filter)
          flt = &(*filter)[0];

        ib = mat_map.begin ();
        ie = mat_map.end ();
        nr = ib->second->get_n_rows ();
        nb = ib->second->get_n_block_size ();
        b_sqr = nb * nb;
        flg.resize (nr);

        if (filter && (t_long)(filter->size ()) != nr)
          {
            // FIXME: 
            if (filter->size () == 0)
              {
                filter->init (nr, 1);
                flt = filter->data ();
              }
            else
              {
                // TODO: FIX
                throw "Filter size problem";
              }
          }
        // create matrix
        m = BS_KERNEL.create_object ("bcsr_matrix");
        m->init (nr, nr, nb, 0);
        rows = &(*(m->get_rows_ptr ()))[0];
        memset (rows, 0, sizeof (t_long) * (nr + 1));

        // step 1: init matrix structure
        for (t_long i = 0; i < nr; ++i)
          flg[i] = -1;

        count = 0;
        for (t_long i = 0; i < nr; ++i)
          {
            // add diagonal into the first place
            flg[i] = count;
            ++count;

            if (flt && flt[i])
              {
                for (it = ib; it != ie; ++it)
                  {
                    i_rows = &(*(it->second->get_rows_ptr ()))[0];
                    i_cols = &(*(it->second->get_cols_ind ()))[0];
                    j1 = i_rows[i];
                    j2 = i_rows[i + 1];
                    for (j = j1; j < j2; ++j)
                      {
                        cl = i_cols[j];
                        if (flt && flt[cl])
                          {
                            // check if element already exist in matrix structure
                            if (flg[cl] < rows[i])
                              {
                                // add new element to the matrix structure
                                flg[cl] = count;
                                ++count;
                              }
                          }
                      }
                  }
              }
            rows[i + 1] = count;
          }
        // step 2: fill matrix
        for (t_long i = 0; i < nr; ++i)
          flg[i] = -1;

        m->alloc_cols_ind_and_values (count, nb);
        cols = &(*(m->get_cols_ind ()))[0];
        values = &(*(m->get_values ()))[0];
        count = 0;
        for (t_long i = 0; i < nr; ++i)
          {
            // add diagonal into the first place
            cols[count] = i;
            flg[i] = count;
            ++count;

            if (flt && flt[i])
              {
                for (it = ib; it != ie; ++it)
                  {
                    i_rows = &(*(it->second->get_rows_ptr ()))[0];
                    i_cols = &(*(it->second->get_cols_ind ()))[0];
                    i_values = &(*(it->second->get_values ()))[0];

                    j1 = i_rows[i];
                    j2 = i_rows[i + 1];
                    for (j = j1; j < j2; ++j)
                      {
                        cl = i_cols[j];
                        if (flt && flt[cl])
                          {
                            // check if element already exist in matrix structure
                            if (flg[cl] < rows[i])
                              {
                                // add new element to the matrix structure
                                cols[count] = cl;
                                flg[cl] = count;
                                memcpy (values + count * b_sqr, i_values + j * b_sqr,
                                        sizeof (t_float) * b_sqr);
                                ++count;
                              }
                            else
                              {
                                // element already exist in matrix
                                MM_SUM (nb, i_values + j * b_sqr, values + flg[cl] * b_sqr);
                              }
                          }
                      }
                  }
              }
            if (rows[i + 1] != count)
              {
                // TODO: FIX
                throw "Internal check error";
              }
          }
      }
    else
      {
        // TODO: FIX
        throw "Internal check error: empty matrix";
      }
    return m;
  }

  BLUE_SKY_TYPE_STD_CREATE (mbcsr_matrix);
  BLUE_SKY_TYPE_STD_COPY (mbcsr_matrix);


  BLUE_SKY_TYPE_IMPL (mbcsr_matrix, mbcsr_matrix_iface, "mbcsr_matrix", "Multi BCSR Matrix class", "Realization of Multi BCSR Matrix");
  }
