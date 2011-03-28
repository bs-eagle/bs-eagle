/** 
 * @file direct_pbuild.cpp
 * @brief 
 * @author 
 * @version 
 * @date 2010-03-16
 */
//#include "amg_stdafx.h"
#include "direct_pbuild.h"

namespace blue_sky
{
  direct_pbuild::direct_pbuild (bs_type_ctor_param) 
                 : amg_pbuild_iface ()
    {
    }
  direct_pbuild::direct_pbuild (const this_t & /*src*/) : bs_refcounter ()
     {
     }
  
  int 
  direct_pbuild::build (sp_bcsr_t matrix, 
                       const t_long n_coarse_size,
                       const t_long /*max_connections*/,
                       spv_long sp_cf_markers,
                       spv_long sp_s_markers,
                       sp_bcsr_t p_matrix)
    {
    t_long n, i, j1, j2, j, cl, j_ind;;
    t_long c_counter;
    t_double sum_d, sum_c, sum_all;
    t_double alpha;
    const double eps = 1.0e-12;

    n = matrix->get_n_rows ();
    t_long *a_cols_ind = &(*(matrix->get_cols_ind ()))[0];
    t_long *a_rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
    t_float *a_values = &(*(matrix->get_values ()))[0];


    if (p_matrix->alloc_rows_ptr (n))
      {
        //TODO: print error message
        BS_ASSERT (false&&"memory allocation failed (bcsr_matrix->alloc_rows_ptr)");
        return -2;
      }

    t_long *p_rows_ptr = &(*(p_matrix->get_rows_ptr ()))[0];
    p_matrix->get_rows_ptr ()->assign (0);

    t_long *cf_markers = &(*sp_cf_markers)[0];
    t_long *s_markers = &(*sp_s_markers)[0];

    c_counter = 0;
    for (i = 0; i < n; ++i)
      {
        if (cf_markers[i] > 0)
          cf_markers[i] = (++c_counter);
      }

    if (c_counter < 1 || n_coarse_size != c_counter)
      return -1;

    // calculate number of elements per row
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
                    if (cf_markers[a_cols_ind[j]] > 0)
                      ++p_rows_ptr[i + 1];
                  }
              }
          }
      }
    // fill p_rows_ptr
    for (i = 2; i <= n; ++i)
      p_rows_ptr[i] += p_rows_ptr[i - 1];

    // allocate memory for P matrix
    //p_matrix->n_block_size = 1;
    p_matrix->set_n_cols (n_coarse_size);
    if (p_matrix->alloc_cols_ind_and_values (p_rows_ptr[n], 1 /*block_size*/))
      {
        // print error message
        return -3;
      }

    t_long *p_cols_ind = &(*(p_matrix->get_cols_ind ()))[0];
    t_float *p_values = &(*(p_matrix->get_values ()))[0];

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
            //sum_d = a_values[a_diag_ind[i]];
            sum_d = a_values[j1];  // diagonal should be the first value in row
            sum_c = 0.0;
            sum_all = 0.0;

            if (sum_d < 0)
              {
                //BOSERR (section::amg, level::error) << "DD ERROR " << i << bs_end;
                //TODO: print error message
                return -1;
              }
            for (j = j1 + 1; j < j2; ++j) // skip diagonal element 
              {
                cl = a_cols_ind[j];
                if (s_markers[j] && cf_markers[cl] > 0)
                  {
                    sum_c += a_values[j];
                  }
                if (a_values[j] < 0)
                  {
                    sum_all += a_values[j];
                  }
                else
                  {
                    sum_d += a_values[j];
                  }
              }
            if (sum_c >= 0 && p_rows_ptr[i + 1] - p_rows_ptr[i] > 0)
              {
                //BOSWARN (section::amg, level::warning) << "SUM C ERROR " << i << "\t " << p_rows_ptr[i + 1] - p_rows_ptr[i] << bs_end;
                //TODO: print error message
              }
            else if (sum_d < eps)
              {
                //BOSWARN (section::amg, level::warning) << "SUM D ERROR " << i << bs_end;
                //TODO: print error message
              }
            else if (p_rows_ptr[i + 1] - p_rows_ptr[i] > 0)
              {
                alpha = -sum_all / sum_c;
                sum_d = 1.0 / sum_d;
                j_ind = p_rows_ptr[i];
                for (j = j1; j < j2; ++j)
                  {
                    cl = a_cols_ind[j];
                    if (s_markers[j] && cf_markers[cl] > 0)
                      {
                        p_values[j_ind] = alpha * a_values[j] * sum_d;
                        p_cols_ind[j_ind] = cf_markers[cl] - 1;
                        ++j_ind;
                      }
                  }
                if (j_ind != p_rows_ptr[i + 1])
                  {
                    //BOSERR (section::amg, level::error) << "ERROR ROWS " << i << bs_end;
                    //TODO: print error message
                    return -1;
                  }
              }

          }

      }
    return 0;
    }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (direct_pbuild);
  BLUE_SKY_TYPE_STD_COPY (direct_pbuild);

  BLUE_SKY_TYPE_IMPL (direct_pbuild, amg_pbuild_iface, "direct_pbuild", "Direct prolangation matrix builder class", "Realization of Direct prolangation matrix builder");
}  // blue_sky namespace
