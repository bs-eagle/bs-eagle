/**
 * @file standart_pbuild.cpp
 * @brief implementation of standart interpolation
 * @author
 * @version
 * @date 2010-03-16
 */
//#include "amg_stdafx.h"
#include "standart_pbuild.h"

namespace blue_sky
{

  standart_pbuild::standart_pbuild (bs_type_ctor_param)
                    : amg_pbuild_iface (),
                      sp_markers (BS_KERNEL.create_object (v_long::bs_type ()))
    {
    }

  standart_pbuild::standart_pbuild (const this_t & /*src*/)
                    : bs_refcounter (),
                      sp_markers (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }

   int
  standart_pbuild::build (sp_bcsr_t matrix,
                                 const t_long n_coarse_size,
                                 const t_long /*max_connections*/,
                                 spv_long sp_cf_markers,
                                 spv_long sp_s_markers,
                                 sp_bcsr_t p_matrix)
    {
    t_long i, n, j, j1, j2, cl;
    t_long jj, jj1, jj2, cl2;
    t_long f_ind;

    t_double sum, diag, d;
    t_long sgn;
    t_long n_threads = 1, thread_num = 0, row_begin, row_end;
    t_long c_counter;
    t_long *markers = 0;

    n = matrix->get_n_rows ();
    t_long *a_cols_ind = &(*(matrix->get_cols_ind ()))[0];
    t_long *a_rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
    t_float *a_values  = &(*(matrix->get_values ()))[0];
    t_long *cf_markers = &(*sp_cf_markers)[0];
    t_long *s_markers  = &(*sp_s_markers)[0];

    if (p_matrix->alloc_rows_ptr (n))
      {
        //TODO: print error message
        BS_ASSERT (false&&"memory allocation failed (bcsr_matrix->alloc_rows_ptr)");
        return -2;
      }

    t_long *p_rows_ptr = &(*(p_matrix->get_rows_ptr ()))[0];
    memset (p_rows_ptr, 0, (n + 1) * sizeof (t_long));

#ifdef BUILD_P_PARALLEL
    n < omp_get_max_threads () ? n_threads = 1 : n_threads = omp_get_max_threads ();
#endif //BUILD_P_PARALLEL

    // local memory allocation
    sp_markers->resize (n * n_threads);
    markers = &(*sp_markers)[0];
    memset (markers, 0, sizeof (t_long) * n * n_threads);

    c_counter = 0;
    for (i = 0; i < n; ++i)
      {
        if (cf_markers[i] > 0)
          cf_markers[i] = (++c_counter);
      }
    if (n_coarse_size < 1 || n_coarse_size != c_counter)
      return -1;
    ///////////////////////
    // First step
    // Calculate P matrix size and allocate memory for it
    ///////////////////////
    // Calculate number of non zeros in row of P matrix
#ifdef BUILD_P_PARALLEL
#pragma omp parallel for private (j, j1, j2, cl)
#endif //BUILD_P_PARALLEL
    for (i = 0; i < n; ++i)
      {
        if (cf_markers[i] > 0) // C point
          ++p_rows_ptr[i + 1]; // only one nonzero value in row

        else // F point
          {
            j1 = a_rows_ptr[i];
            j2 = a_rows_ptr[i + 1];
            // find C points which strongly influens i point
            for (j = j1; j < j2; ++j)
              {
                if (s_markers[j])
                  {
                    cl = a_cols_ind[j];
                    if (cf_markers[cl] > 0)
                      ++p_rows_ptr[i + 1];
                  }
              }
          }
      }
    // fill p_rows_ptr
    for (i = 2; i <= n; ++i)
      p_rows_ptr[i] += p_rows_ptr[i - 1];

    // allocate memory for P matrix
    if (p_matrix->alloc_cols_ind_and_values (p_rows_ptr[n], 1))
      return -3;
    p_matrix->set_n_cols (n_coarse_size);
    t_long *p_cols_ind = &(*(p_matrix->get_cols_ind ()))[0];
    t_float *p_values  = &(*(p_matrix->get_values ()))[0];

    // fill rows for C points and F points with type 0
#ifdef BUILD_P_PARALLEL
#pragma omp parallel private (markers, thread_num, i, j, j1, j2, f_ind, cl, diag, sum, sgn,   \
                              jj1, jj2, jj, d, cl2, row_begin, row_end)
#endif //BUILD_P_PARALLEL
    {
#ifdef BUILD_P_PARALLEL
      thread_num = omp_get_thread_num ();
      markers = &(*sp_markers)[0] + n * thread_num;
#else
      markers = &(*sp_markers)[0];
#endif //BUILD_P_PARALLEL
      if (thread_num < n_threads) // thead is used
        {
          row_begin = n * thread_num / n_threads;
          row_end = n * (thread_num + 1) / n_threads;
        }
      else
        row_begin = row_end = 0;
      for (i = row_begin; i < row_end; ++i)
        {
          if (cf_markers[i] > 0) // C point
            {
              j = p_rows_ptr[i]; // only one row
              p_cols_ind[j] = cf_markers[i] - 1;
              p_values[j] = 1.0;
            }
          else // Interpolation for F points
            {
              f_ind = p_rows_ptr[i];
              j1 = a_rows_ptr[i];
              j2 = a_rows_ptr[i + 1];
              diag = 0.0;
              for (j = j1; j < j2; ++j)
                {
                  if (s_markers[j])
                    {
                      cl = a_cols_ind[j];
                      if (cf_markers[cl] > 0)
                        {
                          markers[cl] = f_ind;
                          p_cols_ind[f_ind] = cf_markers[cl] - 1;
                          p_values[f_ind] = a_values[j];

                          ++f_ind;
                        }
                    }
                  else
                    diag += a_values[j];
                }
              for (j = j1; j < j2; ++j)
                {
                  if (s_markers[j])
                    {
                      cl = a_cols_ind[j];
                      // strong F point
                      if (cf_markers[cl] < 0)
                        {
                          sum = 0.0;
                          sgn = 1;
                          jj1 = a_rows_ptr[cl];
                          jj2 = a_rows_ptr[cl + 1];
                          if (a_values[jj1] < 0)
                            sgn = -1;
                          // for row cl in A matrix
                          for (jj = jj1; jj < jj2; ++jj)
                            {
                              cl2 = a_cols_ind[jj];
                              if (markers[cl2] >= p_rows_ptr[i] && sgn * a_values[jj] < 0)
                                sum += a_values[jj];
                            }
                          if (sum != 0)
                            {
                              d = a_values[j] / sum;
                              for (jj = jj1; jj < jj2; ++jj)
                                {
                                  cl2 = a_cols_ind[jj];
                                  if (markers[cl2] >= p_rows_ptr[i] && sgn * a_values[jj] < 0)
                                    p_values[markers[cl2]] += d * a_values[jj];
                                }
                            }
                          else
                            diag += a_values[j];
                        }
                    }
                }
              diag = 1.0 / diag;
              j1 = p_rows_ptr[i];
              j2 = p_rows_ptr[i + 1];

              for (j = j1; j < j2; ++j)
                {
                  p_values[j] *= -diag;
                  //if (p_values[j] < 0)
                  //  neg_c++;
                  //else if (p_values[j] > 1)
                  //  pos_c++;
                }

            }
          // debug
          //j1 = p_rows_ptr[i];
          //j2 = p_rows_ptr[i + 1];
          //for (j = j1; j < j2; ++j)
          //  printf ("%d %d %d %30.20lf\n", n, i, p_cols_ind[j], p_values[j]);
        }
    } //end parallel (if any)
    // fill rows for F points

    //printf ("NEG %d\t POS %d\n", neg_c, pos_c);


    //compress_p (0.2);
    return 0;
  }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (standart_pbuild);
  BLUE_SKY_TYPE_STD_COPY (standart_pbuild);

  BLUE_SKY_TYPE_IMPL (standart_pbuild, amg_pbuild_iface, "standart_pbuild", "Standart prolongation matrix builder class", "Realization of Standart prolongation matrix builder");
}  // blue_sky namespace
