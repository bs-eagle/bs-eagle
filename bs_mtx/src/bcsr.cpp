/**
 * @file bcsr.cpp
 * @brief
 * @author Oleg Borschuk
 * @date 2009-08-18
 */

#include "bcsr.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  bcsr::bcsr (bs_type_ctor_param)
        : bcsr_amg_matrix_iface (),
        values (BS_KERNEL.create_object (v_float::bs_type ())),
        cols_ind (BS_KERNEL.create_object (v_long::bs_type ())),
        rows_ptr (BS_KERNEL.create_object (v_long::bs_type ()))

    {
      n_rows = n_cols = n_block_size = 0;
    }
  bcsr::bcsr (const bcsr & /*src*/) : bs_refcounter (),
        values (BS_KERNEL.create_object (v_float::bs_type ())),
        cols_ind (BS_KERNEL.create_object (v_long::bs_type ())),
        rows_ptr (BS_KERNEL.create_object (v_long::bs_type ()))
     {
      n_rows = n_cols = n_block_size = 0;
     }


  int
  bcsr::matrix_vector_product_t (spv_double v_, spv_double r_) const
  {
    t_long i, j1, j2, j, cl;

    t_long *r_ptr = &(*rows_ptr)[0];
    t_long *c_ind = &(*cols_ind)[0];
    t_float *val = &(*values)[0];
    t_double *v = &(*v_)[0];
    t_double *r = &(*r_)[0];

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
        j1 = r_ptr[i];
        j2 = r_ptr[i + 1];

        ////////////////////////////////////////////
        // loop through rows in transform matrix
        ////////////////////////////////////////////
        for (j = j1; j < j2; ++j)
          {
            cl = c_ind[j];
            r[cl] += v[i] * val[j];
          }
      }
    return 0;
  }

  // calculate linear combination with current matrix A:
  // r = alpha * A * u + beta * v
  // 1 step: r = beta * v / alpha
  // 2 step: r += alpha * matrix * u
  int
  bcsr::calc_lin_comb (t_double alpha,
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
        //t_long n2 = n - (n % 4);
        // TODO:
        for (; i < n; ++i)
          {
            r[i] = v[i] * d;
          }
      }
    else
      {
        memset (r, 0, sizeof (t_double) * n);
        //r.assign (n, 0);
        //assign (r, 0);
      }

    if (fabs (alpha) > eps)
      {
        matrix_vector_product_internal <true> (u_, r_, alpha);
      }

    return r_code;
  }

  t_double
  bcsr::get_allocated_memory_in_mbytes () const
    {
      t_double d = 0;

      d += sizeof (this);
      d += sizeof (t_float) * values->size ();
      d += sizeof (t_long) * rows_ptr->size ();
      d += sizeof (t_long) * cols_ind->size ();
      d /= 1024 * 1024;
      return d;
    }

  int
  bcsr::init (const t_long new_n_rows,
              const t_long new_n_cols,
              const t_long new_n_block_size,
              const t_long new_n_non_zeros)
  {
    if (new_n_block_size < 1 || new_n_cols < 1)
      return -1;

    n_block_size = new_n_block_size;
    n_cols = new_n_cols;

    if (alloc_rows_ptr (new_n_rows)
        || alloc_cols_ind_and_values (new_n_non_zeros, new_n_block_size))
      {
        return -2;
      }

    return 0;
  }

  int
  bcsr::alloc_rows_ptr (const t_long new_n_rows)
  {
    if (new_n_rows < 1)
      return -1;

    n_rows = new_n_rows;

    rows_ptr->resize (n_rows + 1);
    memset (&(*rows_ptr)[0], 0, sizeof (t_long) * (n_rows + 1));

    return 0;
  }

  int
  bcsr::alloc_cols_ind (const t_long new_n_non_zeros)
  {
    t_long nnz = new_n_non_zeros;

    if (nnz < 1)
      nnz = 1;

    cols_ind->resize (nnz);
    memset (&(*cols_ind)[0], -1, sizeof (t_long) * nnz);
    return 0;
  }

  int
  bcsr::alloc_values (const t_long new_n_non_zeros,
                               const t_long new_n_block_size)
  {
    t_long b_sqr = new_n_block_size * new_n_block_size;
    t_long nnz = new_n_non_zeros;

    n_block_size = new_n_block_size;
    if (nnz < 1)
      nnz = 1;
    if (b_sqr < 1)
      b_sqr = 1;

    values->resize (b_sqr * nnz);
    memset (&(*values)[0], 0, sizeof (t_float) * nnz * b_sqr);

    return 0;
  }

  int
  bcsr::alloc_cols_ind_and_values (const t_long new_n_non_zeros,
                                   const t_long new_n_blok_size)
  {
    if (alloc_cols_ind (new_n_non_zeros) || alloc_values (new_n_non_zeros, new_n_blok_size))
      return -2;

    return 0;
  }

  int
  bcsr::init_by_matrix  (sp_bcsr_matrix_iface_t m)
    {
      rows_ptr = m->get_rows_ptr ()->clone ();
      cols_ind = m->get_cols_ind ()->clone ();
      values = m->get_values ()->clone ();
      n_rows = m->get_n_rows ();
      n_cols = m->get_n_cols ();
      n_block_size = m->get_n_block_size ();
      return 0;
    }

  int
  bcsr::init_struct (const t_long new_n_rows,
                     const t_long new_n_cols,
                     const t_long new_n_non_zeros)
{
  if (new_n_cols < 1)
    return -1;

  n_cols = new_n_cols;
  if (alloc_rows_ptr (new_n_rows) || alloc_cols_ind (new_n_non_zeros))
    return -2;

  return 0;
}

  int
  bcsr::build_transpose_struct (sp_bcsr_matrix_iface_t m,
                                const t_long rows_offset,
                                const t_long cols_offset,
                                const t_long new_n_rows)
  {

    t_long i,j;
    t_long row_ind, j1, j2;
    t_long nnz;
    t_long nr;
    t_long *r_ptr;
    t_long *c_ind;
    t_long *in_rows_ptr = &(*(m->get_rows_ptr ()))[0];
    t_long *in_cols_ind = &(*(m->get_cols_ind ()))[0];
    //t_float *v_ptr;

    nnz = m->get_n_non_zeros ();
    t_long in_n_rows = m->get_n_rows ();
    t_long in_n_cols = m->get_n_cols ();
    t_long in_n_block_size = m->get_n_block_size ();

    if (new_n_rows == 0)
      nr = in_n_cols;
    else
      nr = new_n_rows;

    if (init (nr, in_n_rows, in_n_block_size, nnz))
      return -1;

    r_ptr = &(*rows_ptr)[0];
    c_ind = &(*cols_ind)[0];
    // Count the number of entries in each column of matrix (row of transposed matrix) and fill the rows_ptr array.

    for (i = 0; i < nnz; i++)
      {
        ++r_ptr[in_cols_ind[i] + 1 - rows_offset];
      }

    for (i = 2; i <= n_rows; i++)
      {
        r_ptr[i] += r_ptr[i - 1];
      }

    // Load the values and column numbers of transposed matrix

    for (i = 0; i < in_n_rows; i++)
      {
        j1 = in_rows_ptr[i];
        j2 = in_rows_ptr[i + 1];
        for (j = j1; j < j2; j++)
          {
            row_ind = in_cols_ind[j] - rows_offset;
            c_ind[r_ptr[row_ind]] = i + cols_offset;
            r_ptr[row_ind]++;
          }
      }
    // rows_ptr now points to the *end* of the jth row of entries
    // instead of the beginning.  Restore rows_ptr to front of row.

    for (i = n_rows; i > 0; i--)
      {
        r_ptr[i] = r_ptr[i - 1];
      }

    r_ptr[0] = 0;

    return 0;
  }

  int
  bcsr::build_transpose (sp_bcsr_matrix_iface_t m,
                         const t_long rows_offset,
                         const t_long cols_offset,
                         const t_long new_n_rows)
  {

    t_long i,j;
    t_long row_ind, j1, j2;
    t_long nnz;
    t_long nr;
    t_long *r_ptr;
    t_long *c_ind;
    t_long *in_rows_ptr = &(*(m->get_rows_ptr ()))[0];
    t_long *in_cols_ind = &(*(m->get_cols_ind ()))[0];
    t_float *v_ptr;
    t_float *in_values = &(*(m->get_values ()))[0];

    nnz = m->get_n_non_zeros ();
    t_long in_n_rows = m->get_n_rows ();
    t_long in_n_cols = m->get_n_cols ();
    t_long in_n_block_size = m->get_n_block_size ();

    if (new_n_rows == 0)
      nr = in_n_cols;
    else
      nr = new_n_rows;

    if (init (nr, in_n_rows, in_n_block_size, nnz))
      return -1;

    r_ptr = &(*rows_ptr)[0];
    c_ind = &(*cols_ind)[0];
    v_ptr = &(*values)[0];
    // Count the number of entries in each column of matrix (row of transposed matrix) and fill the rows_ptr array.

    for (i = 0; i < nnz; i++)
      {
        ++r_ptr[in_cols_ind[i] + 1 - rows_offset];
      }

    for (i = 2; i <= n_rows; i++)
      {
        r_ptr[i] += r_ptr[i - 1];
      }

    // Load the values and column numbers of transposed matrix

    for (i = 0; i < in_n_rows; i++)
      {
        j1 = in_rows_ptr[i];
        j2 = in_rows_ptr[i + 1];
        for (j = j1; j < j2; j++)
          {
            row_ind = in_cols_ind[j] - rows_offset;
            c_ind[r_ptr[row_ind]] = i + cols_offset;
            v_ptr[r_ptr[row_ind]] = in_values[j];
            r_ptr[row_ind]++;
          }
      }
    // rows_ptr now points to the *end* of the jth row of entries
    // instead of the beginning.  Restore rows_ptr to front of row.

    for (i = n_rows; i > 0; i--)
      {
        r_ptr[i] = r_ptr[i - 1];
      }

    r_ptr[0] = 0;

    return 0;
  }
  int
  bcsr::triple_matrix_product (sp_bcsr_matrix_iface_t r_matrix, sp_bcsr_matrix_iface_t a_matrix,
                               sp_bcsr_matrix_iface_t p_matrix, const bool update)
  {
    t_long i, i1, i2, j;
    t_long jr1, jr2, jr;
    t_long ja1, ja2, ja;
    t_long jp1, jp2, jp;
    t_long index_counter, i_row_ptr; // indexes of RAP values - current element and current row start element respectively
    t_long thread_num = 0, n_threads = 1, row_begin, row_end, row_begin_s, row_end_s;
    t_float r_a_product;
    t_long *r_ptr;
    t_long *c_ind;
    t_float *v_ptr;

    t_long *r_rows_ptr = &(*(r_matrix->get_rows_ptr ()))[0];
    t_long *r_cols_ind = &(*(r_matrix->get_cols_ind ()))[0];
    t_float *r_values   = &(*(r_matrix->get_values ()))[0];
    t_long  r_n_rows   = r_matrix->get_n_rows ();

    t_long *a_rows_ptr = &(*(a_matrix->get_rows_ptr ()))[0];
    t_long *a_cols_ind = &(*(a_matrix->get_cols_ind ()))[0];
    t_float *a_values   = &(*(a_matrix->get_values ()))[0];
    t_long  a_n_cols   = a_matrix->get_n_cols ();

    t_long *p_rows_ptr = &(*(p_matrix->get_rows_ptr ()))[0];
    t_long *p_cols_ind = &(*(p_matrix->get_cols_ind ()))[0];
    t_float *p_values   = &(*(p_matrix->get_values ()))[0];
    t_long  p_n_cols   = p_matrix->get_n_cols ();
    t_long block_size  = a_matrix->get_n_block_size ();

    // works only with non-blocked matrices
    if (block_size > 1)
      {
        //TODO: print error message
        return -2;
      }

    if (!update)
      {
        // initializing RAP matrix
        //result->n_block_size = A->n_block_size;
        //result->n_cols = P->n_cols;
        if (alloc_rows_ptr (r_n_rows))
          {
            //TODO: print error message
            return -1;
          }
      }
    r_ptr = &(*rows_ptr)[0];

    // initializing some stuff (c)

#ifdef TRIPLE_MATRIX_PRODUCT_PARALLEL
#pragma omp master
    {
      r_matrix->get_n_rows () < omp_get_max_threads () ? n_threads = r_matrix->get_n_rows () : n_threads = omp_get_max_threads ();
    }
#endif //TRIPLE_MATRIX_PRODUCT_PARALLEL

    a_marker_vec.assign (a_n_cols * n_threads, -1 );
    p_marker_vec.assign (p_n_cols * n_threads, -1 );
    t_long *a_marker = &a_marker_vec[0];
    t_long *p_marker = &p_marker_vec[0];

    r_ptr[0] = 0;

    if (!update)
      {
        ///////////////////////////////////////
        //  STEP 1 - calculating non-zero RAP elements count
        ///////////////////////////////////////
        // loop through R rows
#ifdef TRIPLE_MATRIX_PRODUCT_PARALLEL
#pragma omp parallel  firstprivate (A_marker, P_marker)                                \
                          private (i, thread_num, row_begin, row_end, i_row_ptr, jr1, jr2, jr, ja1, ja2, ja, \
                                   jp1, jp2, jp, i1, i2, j, index_counter)
#endif //TRIPLE_MATRIX_PRODUCT_PARALLEL
        {

#ifdef TRIPLE_MATRIX_PRODUCT_PARALLEL
          thread_num = omp_get_thread_num ();
          a_marker += a_matrix->get_n_cols () * thread_num;
          p_marker += p_matrix->get_n_cols () * thread_num;
#endif //TRIPLE_MATRIX_PRODUCT_PARALLEL
          index_counter = 0;
          if (thread_num < n_threads) // thead is used
            {
              row_begin = r_n_rows * thread_num / n_threads;
              row_end = r_n_rows * (thread_num + 1) / n_threads;
            }
          else
            row_begin = row_end = 0;

          for (i = row_begin; i < row_end; i++)
            {
              i_row_ptr = index_counter;

              jr1 = r_rows_ptr[i];
              jr2 = r_rows_ptr[i + 1];
              // loop through R columns in i row
              for (jr = jr1; jr < jr2; jr++)
                {
                  i1 = r_cols_ind[jr];

                  ja1 = a_rows_ptr[i1];
                  ja2 = a_rows_ptr[i1 + 1];
                  // loop through A columns in i1 row
                  for (ja = ja1; ja < ja2; ja++)
                    {
                      i2 = a_cols_ind[ja];

                      //  Check A_marker to see if point i2 has been previously
                      //  visited. New entries in RAP only occur from unmarked points.

                      if (a_marker[i2] != i)
                        {

                          //  Mark i2 as visited.

                          a_marker[i2] = i;

                          jp1 = p_rows_ptr[i2];
                          jp2 = p_rows_ptr[i2 + 1];
                          // loop through P columns in i2 row
                          for (jp = jp1; jp < jp2; jp++)
                            {
                              j = p_cols_ind[jp];

                              //  Check P_marker to see that RAP{i,j} has not already
                              //  been accounted for. If it has not, mark it and increment
                              //  counter.

                              if (p_marker[j] < i_row_ptr)
                                {
                                  // we have new entry in RAP matrix
                                  p_marker[j] = index_counter;
                                  index_counter++;
                                }
                            }
                        }
                    }
                }
              r_ptr[i + 1] = index_counter;
            }
#ifdef TRIPLE_MATRIX_PRODUCT_PARALLEL
#pragma omp barrier
#pragma omp single
#endif //TRIPLE_MATRIX_PRODUCT_PARALLEL
          // collect results of all threads at the head of each thread's chunk
          {
            i_row_ptr = 0;
            for (j = 1; j < n_threads - 1; j++)
              {
                row_begin_s = r_n_rows * j / n_threads;
                row_end_s = r_n_rows * (j + 1) / n_threads;
                r_ptr[row_end_s] += r_ptr[row_begin_s];
              }
          }
          // each thread adds collected head value to each member of own chunk
          if (thread_num != 0)
            {
              for (i = row_begin; i < row_end - 1; i++)
                {
                  r_ptr[i + 1] += r_ptr[row_begin];
                }
              if (row_end == r_n_rows)
                r_ptr[row_end] += r_ptr[row_begin];
            }
        } // end omp parallel



        // now we have (RAP_rows_ptr[result->n_rows]) non-zero values in RAP matrix.
        // Let`s finish it`s initialization
        if (alloc_cols_ind_and_values(r_ptr[n_rows], 1))
          //TODO: print error message
          return -1;
        set_n_cols (n_rows);
      }
    c_ind = &(*cols_ind)[0];
    v_ptr = &(*values)[0];

    // initializing some stuff (c)


    a_marker_vec.assign (a_n_cols * n_threads, -1 );
    p_marker_vec.assign (p_n_cols * n_threads, -1 );
    a_marker = &a_marker_vec[0];
    p_marker = &p_marker_vec[0];



    ///////////////////////////////////////
    //  STEP 2 - calculating a product
    ///////////////////////////////////////
#ifdef TRIPLE_MATRIX_PRODUCT_PARALLEL
#pragma omp parallel  firstprivate (a_marker, p_marker)                                   \
                      private (i, thread_num, row_begin, row_end, i_row_ptr, jr1, jr2, jr, ja1, ja2, ja, \
                               jp1, jp2, jp, i1, i2, j, index_counter, r_a_product)
#endif //TRIPLE_MATRIX_PRODUCT_PARALLEL
    {
#ifdef TRIPLE_MATRIX_PRODUCT_PARALLEL
      thread_num = omp_get_thread_num ();
      a_marker += a_matrix.get_n_cols () * thread_num;
      p_marker += p_matrix.get_n_cols () * thread_num;
#endif //TRIPLE_MATRIX_PRODUCT_PARALLEL

      if (thread_num < n_threads) // thead is used
        {
          row_begin = r_n_rows * thread_num / n_threads;
          row_end = r_n_rows * (thread_num + 1) / n_threads;
        }
      else
        row_begin = row_end = 0;
      index_counter = r_ptr[row_begin];
      for (i = row_begin; i < row_end; i++)
        {
          i_row_ptr = index_counter;
          p_marker[i] = index_counter;
          v_ptr[index_counter] = 0.0;
          c_ind[index_counter] = i;
          ++index_counter;


          jr1 = r_rows_ptr[i];
          jr2 = r_rows_ptr[i + 1];

          // loop through R columns in i row
          for (jr = jr1; jr < jr2; jr++)
            {
              i1 = r_cols_ind[jr];

              ja1 = a_rows_ptr[i1];
              ja2 = a_rows_ptr[i1 + 1];

              // loop through A columns in i1 row
              for (ja = ja1; ja < ja2; ja++)
                {
                  i2 = a_cols_ind[ja];

                  r_a_product = r_values[jr] * a_values[ja];

                  //  Check A_marker to see if point i2 has been previously
                  //  visited. New entries in RAP only occur from unmarked points.

                  if (a_marker[i2] != i)
                    {

                      //  Mark i2 as visited.

                      a_marker[i2] = i;

                      jp1 = p_rows_ptr[i2];
                      jp2 = p_rows_ptr[i2 + 1];

                      for (jp = jp1; jp < jp2; jp++)
                        {
                          j = p_cols_ind[jp];

                          //  Check P_marker to see that RAP{i,j} has not already
                          //  been accounted for. If it has not, mark it and increment
                          //  counter.

                          if (p_marker[j] < i_row_ptr)
                            {
                              // we have new entry in RAP matrix
                              v_ptr[index_counter] = p_values[jp] * r_a_product;
                              c_ind[index_counter] = j;

                              p_marker[j] = index_counter;
                              index_counter++;
                            }
                          else
                            {
                              // we have already visited this element, upgrade it
                              v_ptr[p_marker[j]] += p_values[jp] * r_a_product;
                            }
                        }
                    }
                  else
                    {
                      //  i2 is already visited. No new RAP entries. Just add contribution.

                      jp1 = p_rows_ptr[i2];
                      jp2 = p_rows_ptr[i2 + 1];

                      for (jp = jp1; jp < jp2; jp++)
                        {
                          j = p_cols_ind[jp];
                          v_ptr[p_marker[j]] += p_values[jp] * r_a_product;
                        }
                    }
                }
            }
        }
    } // end parallel

    return 0;
  }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (bcsr);
  BLUE_SKY_TYPE_STD_COPY (bcsr);

  BLUE_SKY_TYPE_IMPL (bcsr, bcsr_amg_matrix_iface, "bcsr_matrix", "Block CSR Matrix class", "Realization of Block CSR Matricies");
}  // blue_sky namespace
