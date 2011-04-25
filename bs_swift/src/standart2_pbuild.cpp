/**
 * @file standart2_pbuild.cpp
 * @brief Stuben's standart interpolation (extended interpolation, using F-F connections)
 * @author
 * @version
 * @date 2010-03-16
 */
//#include "amg_stdafx.h"
#include "standart2_pbuild.h"

namespace blue_sky
{

  standart2_pbuild::standart2_pbuild (bs_type_ctor_param)
                    : amg_pbuild_iface (),
                      sp_markers (BS_KERNEL.create_object (v_long::bs_type ()))
    {
    }

  standart2_pbuild::standart2_pbuild (const this_t & /*src*/)
                    : bs_refcounter (),
                      sp_markers (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }

   int
  standart2_pbuild::build (sp_bcsr_t matrix,
                           const t_long n_coarse_size,
                           const t_long /*max_connections*/,
                           spv_long sp_cf_markers,
                           spv_long sp_s_markers,
                           sp_bcsr_t p_matrix)
    {
    t_long n, i, j1, j2, j, cl, j_ind, j_start, jj, jj1, jj2, cl2;
    t_long c_counter;
    t_double sum_d, sum_c, sum_all;
    t_double alpha;

    n = matrix->get_n_rows ();
    t_long  *a_cols_ind = &(*(matrix->get_cols_ind ()))[0];
    t_long  *a_rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
    t_float *a_values   = &(*(matrix->get_values ()))[0];

    if (p_matrix->alloc_rows_ptr (n))
      {
        //TODO: print error message
        assert(false&&"memory allocation failed (bcsr_matrix->alloc_rows_ptr)");
        return -2;
      }

    t_long *p_rows_ptr = &(*(p_matrix->get_rows_ptr ()))[0];
    memset (p_rows_ptr, 0, sizeof (t_long) * (n + 1));

    // local memory allocation
    sp_markers->resize (n);
    sp_markers->assign (-1);

    t_long *markers = &(*sp_markers)[0];
    t_long *cf_markers = &(*sp_cf_markers)[0];
    t_long *s_markers = &(*sp_s_markers)[0];

    c_counter = 0;
    for (i = 0; i < n; ++i)
      {
        if (cf_markers[i] > 0)
          cf_markers[i] = (++c_counter);
      }

    if (n_coarse_size < 1 || n_coarse_size != c_counter)
      return -1;

    // calculate number of elements per row
    j_ind = 0;
    for (i = 0; i < n; ++i)
      {
        if (cf_markers[i] > 0) // C point
          ++p_rows_ptr[i + 1]; // only one nonzero value in row

        else // F point
          {
            j1 = a_rows_ptr[i];
            j2 = a_rows_ptr[i + 1];
            j_start = j_ind;

            // find C points which strongly influens i point
            for (j = j1; j < j2; ++j)
              {
                if (s_markers[j])
                  {
                    cl = a_cols_ind[j];
                    if (cf_markers[cl] > 0)
                      {
                        if (markers[cl] < j_start)
                          {
                            markers[cl] = j_ind;
                            ++p_rows_ptr[i + 1];
                            ++j_ind;
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
                                if (cf_markers[cl2] > 0 && markers[cl2] < j_start)
                                  {
                                    markers[cl2] = j_ind;
                                    ++p_rows_ptr[i + 1];
                                    ++j_ind;
                                  }
                              }
                          }
                      }
                  }
              }
          }
      }
    // fill p_rows_ptr
    for (i = 2; i <= n; ++i)
      p_rows_ptr[i] += p_rows_ptr[i - 1];

    // allocate memory for P matrix
    if (p_matrix->alloc_cols_ind_and_values (p_rows_ptr[n], 1))
      {
        // print error message
        return -3;
      }
    p_matrix->set_n_cols (n_coarse_size);

    t_long *p_cols_ind = &(*(p_matrix->get_cols_ind ()))[0];
    t_float *p_values = &(*(p_matrix->get_values ()))[0];

    //memset (p_values, 0, sizeof (double) * p_rows_ptr[n]);
    sp_markers->assign (-1);

    // fill P matrix
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


            // STAGE 1 fill values from strong F-connections
            for (j = j1; j < j2; ++j)
              {
                cl = a_cols_ind[j];
                if (s_markers[j])
                  {
                    if (cf_markers[cl] > 0) // C-point
                      {
                        if (markers[cl] >= j_start)
                          {
                            p_values[markers[cl]] += a_values[j];
                          }
                        else
                          {
                            markers[cl] = j_ind;
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
                                    if (markers[cl2] >= j_start)
                                      {
                                        p_values[markers[cl2]] += alpha * a_values[jj]
                                                                  * sum_d * a_values[j];
                                      }
                                    else
                                      {
                                        markers[cl2] = j_ind;
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
                // TODO: print error message
                return -1;
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
    //compress_p (0.5);
    return 0;
  }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (standart2_pbuild);
  BLUE_SKY_TYPE_STD_COPY (standart2_pbuild);

  BLUE_SKY_TYPE_IMPL (standart2_pbuild, amg_pbuild_iface, "standart2_pbuild", "standart2 prolongation matrix builder class", "Realization of standart2 prolongation matrix builder");
}  // blue_sky namespace
