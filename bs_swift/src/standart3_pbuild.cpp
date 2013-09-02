/**
 * @file standart3_pbuild.cpp
 * @brief parallel implementation of Stuben's standart interpolation
 * @author Mark Khait
 * @version
 * @date 2010-03-16
 */
//#include "amg_stdafx.h"
#include "standart3_pbuild.h"

namespace blue_sky
{

  standart3_pbuild::standart3_pbuild (bs_type_ctor_param)
                    : amg_pbuild_iface (),
                      sp_cur_connections (BS_KERNEL.create_object (v_long::bs_type ())),
                      sp_cur_connections_markers (BS_KERNEL.create_object (v_long::bs_type ()))

    {
    }

  standart3_pbuild::standart3_pbuild (const this_t & /*src*/)
                    : bs_refcounter (),
                      sp_cur_connections (BS_KERNEL.create_object (v_long::bs_type ())),
                      sp_cur_connections_markers (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }

   int
  standart3_pbuild::build (sp_bcsr_t matrix,
                           const t_long n_coarse_size,
                           const t_long max_connections,
                           spv_long sp_cf_markers,
                           spv_long sp_s_markers,
                           sp_bcsr_t p_matrix)
    {
    t_float *p_values;
    t_long *p_cols_ind;
    t_long n, i, j1, j2, j, cl, j_ind, j_start, jj, jj1, jj2, cl2, k;
    t_long c_counter;
    t_long match;
    t_double sum_d, sum_c, sum_all;
    t_double alpha;
    t_long n_cur_connections;

    // internal checks
    n = matrix->get_n_rows ();
    t_long *a_cols_ind = &(*(matrix->get_cols_ind ()))[0];
    t_long *a_rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
    t_float *a_values  = &(*(matrix->get_values ()))[0];
    t_long *cf_markers = &(*sp_cf_markers)[0];
    t_long *s_markers  = &(*sp_s_markers)[0];
    t_long *cur_connections_markers;
    t_long *cur_connections;

    if (p_matrix->alloc_rows_ptr (n))
      {
        //TODO: print error message
        assert(false&&"memory allocation failed (bcsr_matrix->alloc_rows_ptr)");
        return -2;
      }

    t_long *p_rows_ptr = &(*(p_matrix->get_rows_ptr ()))[0];
    memset (p_rows_ptr, 0, sizeof (t_long) * (n + 1));

    c_counter = 0;
    for (i = 0; i < n; ++i)
      {
        if (cf_markers[i] > 0)
          cf_markers[i] = (++c_counter);
      }
    if (n_coarse_size < 1 || n_coarse_size != c_counter)
      return -1;

// calculate number of elements per row
#ifdef BUILD_P_STANDART_PARALLEL
#pragma omp parallel private (cur_connections, cur_connections_markers, n_cur_connections,i, j, j1, j2, jj, jj1, jj2, \
                              cl, cl2, k, match, j_ind, sum_d, sum_c, sum_all, alpha)
    {
#endif //BUILD_P_STANDART_PARALLEL

      sp_cur_connections->resize (max_connections * max_connections);
      sp_cur_connections_markers->resize (max_connections * max_connections);
      cur_connections_markers = &(*sp_cur_connections_markers)[0];
      cur_connections = &(*sp_cur_connections)[0];
      memset (cur_connections, 0, max_connections * max_connections * sizeof (t_long));
      memset (cur_connections_markers, 0, max_connections * max_connections * sizeof (t_long));

#ifdef BUILD_P_STANDART_PARALLEL
#pragma omp for
#endif //BUILD_P_STANDART_PARALLEL
      for (i = 0; i < n; ++i)
        {
          if (cf_markers[i] > 0) // C point
            ++p_rows_ptr[i + 1]; // only one nonzero value in row

          else // F point
            {
              j1 = a_rows_ptr[i];
              j2 = a_rows_ptr[i + 1];
              n_cur_connections = 0;

              // find C points which strongly influens i point
              for (j = j1; j < j2; ++j)
                {
                  if (s_markers[j])
                    {
                      cl = a_cols_ind[j];
                      if (cf_markers[cl] > 0)
                        {
                          for (k = 0, match = 0; k < n_cur_connections; k++)
                            {
                              if (cur_connections[k] == cl)
                                {
                                  match = 1;
                                  break;
                                }
                            }
                          if (!match)
                            {
                              cur_connections[n_cur_connections] = cl;
                              ++p_rows_ptr[i + 1];
                              ++n_cur_connections;
                            }
                        }
                      else // strong F-point connecton
                        {
                          jj1 = a_rows_ptr[cl];
                          jj2 = a_rows_ptr[cl + 1];
                          for (jj = jj1; jj < jj2; ++jj)
                            {
                              if (s_markers[jj])
                                {
                                  cl2 = a_cols_ind[jj];
                                  if (cf_markers[cl2] > 0)
                                    {
                                      for (k = 0, match = 0; k < n_cur_connections; k++)
                                        {
                                          if (cur_connections[k] == cl2)
                                            {
                                              match = 1;
                                              break;
                                            }
                                        }
                                      if (!match)
                                        {
                                          cur_connections[n_cur_connections] = cl2;
                                          ++p_rows_ptr[i + 1];
                                          ++n_cur_connections;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

#ifdef BUILD_P_STANDART_PARALLEL
#pragma omp single
#endif //BUILD_P_STANDART_PARALLEL
      {
        // fill p_rows_ptr
        for (i = 2; i <= n; ++i)
          p_rows_ptr[i] += p_rows_ptr[i - 1];
      }

#ifdef BUILD_P_STANDART_PARALLEL
#pragma omp single
#endif //BUILD_P_STANDART_PARALLEL
      {
        // allocate memory for P matrix
        if (p_matrix->alloc_cols_ind_and_values (p_rows_ptr[n], 1))
          {
#ifdef BUILD_P_STANDART_PARALLEL
            BOSERR (section::amg, level::error) << "Build p memory allocations error" << bs_end;
#else //BUILD_P_STANDART_PARALLEL
            // print error message
            return -3;
#endif //BUILD_P_STANDART_PARALLEL
          }
        p_matrix->set_n_cols (n_coarse_size);
        p_cols_ind = &(*(p_matrix->get_cols_ind ()))[0];
        p_values = &(*(p_matrix->get_values ()))[0];
      }

      // fill P matrix

#ifdef BUILD_P_STANDART_PARALLEL
#pragma omp for
#endif //BUILD_P_STANDART_PARALLEL
      for (i = 0; i < n; ++i)
        {
          if (cf_markers[i] > 0) // C point
            {
              j = p_rows_ptr[i]; // only one row
              p_cols_ind[j] = cf_markers[i] - 1;
              p_values[j] = 1.0;
            }
          else // interpolation for F points
            {
              j1 = a_rows_ptr[i];
              j2 = a_rows_ptr[i + 1];
              j_ind = j_start = p_rows_ptr[i];
              n_cur_connections = 0;

              // STAGE 1 fill values from strong F-connections
              for (j = j1; j < j2; ++j)
                {
                  cl = a_cols_ind[j];
                  if (s_markers[j])
                    {
                      if (cf_markers[cl] > 0) // C-point
                        {

                          for (k = 0, match = 0; k < n_cur_connections; k++)
                            {
                              if (cur_connections[k] == cl)
                                {
                                  match = 1;
                                  p_values[cur_connections_markers[k]] += a_values[j];
                                  break;
                                }
                            }
                          if (!match)
                            {
                              cur_connections[n_cur_connections] = cl;
                              cur_connections_markers[n_cur_connections] = j_ind;
                              n_cur_connections++;
                              p_values[j_ind] = a_values[j];
                              p_cols_ind[j_ind] = cf_markers[cl] - 1;
                              ++j_ind;
                            }
                        }
                      else // F-point
                        {
                          jj1 = a_rows_ptr[cl];
                          jj2 = a_rows_ptr[cl + 1];
#if 0
                          alpha = -1.0 / a_values[a_diag_ind[cl]];
                          for (jj = jj1; jj < jj2; ++jj)
                            {
                              cl2 = a_cols_ind[jj];
                              if (s_markers[jj] && cf_markers[cl2] > 0)
                                {
                                  if (markers[cl2] >= j_start)
                                    {
                                      p_values[markers[cl2]] += alpha * a_values[jj]
                                                                * a_values[j];
                                    }
                                  else
                                    {
                                      markers[cl2] = j_ind;
                                      p_values[j_ind] = alpha * a_values[jj] * a_values[j];
                                      p_cols_ind[j_ind] = cf_markers[cl2] - 1;
                                      ++j_ind;
                                    }
                                }
                            }
#endif //0

#if 1
                          sum_d = a_values[jj1];
                          sum_c = 0.0;
                          sum_all = 0.0;
                          for (jj = jj1; jj < jj2; ++jj)
                            {
                              cl2 = a_cols_ind[jj];
                              if (cl2 != cl)
                                {
                                  if (s_markers[jj] && cf_markers[cl2] > 0)
                                    {
                                      sum_c += a_values[jj];
                                    }
                                  if (a_values[jj] < 0 && s_markers[jj])
                                    {
                                      sum_all += a_values[jj];
                                    }
                                  else
                                    {
                                      sum_d += a_values[jj];
                                    }
                                }
                            }
                          if (sum_c < 0)
                            {
                              alpha = -sum_all / sum_c;
                              sum_d = 1.0 / sum_d;
                              for (jj = jj1; jj < jj2; ++jj)
                                {
                                  cl2 = a_cols_ind[jj];
                                  if (s_markers[jj] && cf_markers[cl2] > 0)
                                    {
                                      for (k = 0, match = 0; k < n_cur_connections; k++)
                                        {
                                          if (cur_connections[k] == cl2)
                                            {
                                              match = 1;
                                              p_values[cur_connections_markers[k]] += alpha * a_values[jj]
                                                                                      * sum_d * a_values[j];
                                              break;
                                            }
                                        }
                                      if (!match)
                                        {
                                          cur_connections[n_cur_connections] = cl2;
                                          cur_connections_markers[n_cur_connections] = j_ind;
                                          n_cur_connections++;
                                          p_values[j_ind] = alpha * a_values[jj] * sum_d * a_values[j];
                                          p_cols_ind[j_ind] = cf_markers[cl2] - 1;
                                          ++j_ind;
                                        }
                                    }
                                }

                            }
#endif // 0
                        }

                    }
                }
              if (j_ind != p_rows_ptr[i + 1])
                {
                  //BOSERR (section::amg, level::error) << "ERROR ROWS BLA " << i << bs_end;
                  //TODO: print error message
#ifndef BUILD_P_STANDART_PARALLEL
                  return -1;
#endif //BUILD_P_STANDART_PARALLEL
                }

              //sum_c = 0.0;

              jj1 = p_rows_ptr[i];
              jj2 = p_rows_ptr[i + 1];
              sum_d = a_values[a_rows_ptr[i]];//diagonal element
              //for (jj = jj1; jj < jj2; ++jj)
              //  sum_c += p_values[jj];
              //sum_all = sum_c;
              for (j = j1; j < j2; ++j)
                {

                  cl = a_cols_ind[j];
                  if (cl != i)
                    {
                      // if (a_values[j] < 0 && markers[cl] < j_start)
                      //   {
                      //     sum_all += a_values[j];
                      //   }
                      // else if (markers[cl] < j_start)
                      if (s_markers[j] == 0)
                        {
                          sum_d += a_values[j];
                        }
                    }
                }
              if (p_rows_ptr[i + 1] - p_rows_ptr[i] > 0)
                {
                  //alpha = -sum_all / sum_c;
                  sum_d = -1.0 / sum_d;
                  for (jj = jj1; jj < jj2; ++jj)
                    {
                      p_values[jj] *= sum_d;
                    }
                }

            }

        }
#ifdef BUILD_P_STANDART_PARALLEL
    }
#endif //BUILD_P_STANDART_PARALLEL

    //compress_p (0.5);
    return 0;
  }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (standart3_pbuild);
  BLUE_SKY_TYPE_STD_COPY (standart3_pbuild);

  BLUE_SKY_TYPE_IMPL (standart3_pbuild, amg_pbuild_iface, "standart3_pbuild", "standart3 prolongation matrix builder class", "Realization of standart3 prolongation matrix builder");
}  // blue_sky namespace
