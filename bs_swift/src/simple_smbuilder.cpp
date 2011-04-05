/** 
 * @file simple_smbuilder.cpp
 * @brief 
 * @author 
 * @version 
 * @date 2010-03-09
 */
#include "simple_smbuilder.h"

namespace blue_sky
{
  
  simple_smbuilder::simple_smbuilder (bs_type_ctor_param) 
                        : amg_smbuilder_iface ()
    {
    }
  
  simple_smbuilder::simple_smbuilder (const this_t & /*src*/) : bs_refcounter ()
     {
     }
  
   int 
  simple_smbuilder::build (sp_bcsr_t a_matrix, 
                                   t_double strength_threshold, 
                                   t_double max_row_sum,
                                   spv_long sp_s_markers) const
    {
    t_long i, j1, j2, j;
    t_long n, nnz;
    t_double row_max;
    t_double row_sum;
    t_double diag;
    t_long cur_max;
#ifdef BUILD_STRENGTH_MATRIX_PARALLEL
    t_long thread_num;
    i_vector_type_t max_connections_local = max_connections_thread;
    int r_code = 0;
#endif //BUILD_STRENGTH_MATRIX_PARALLEL

    // internal checks and initialization
    assert (a_matrix);

    const t_float *data = &(*(a_matrix->get_values ()))[0];
    const t_long  *rows = &(*(a_matrix->get_rows_ptr ()))[0];

    nnz = a_matrix->get_n_non_zeros ();
    sp_s_markers->resize (nnz);
    sp_s_markers->assign (0);
    t_long *s_markers = &(*sp_s_markers)[0];

    n = a_matrix->get_n_rows ();
    t_long max_connections = 0;
#ifdef BUILD_STRENGTH_MATRIX_PARALLEL
#pragma omp parallel firstprivate (max_connections_local) private (thread_num)
    {
      thread_num = omp_get_thread_num ();
      max_connections_local += thread_num;
      int thread_max = 0;
#pragma omp for private (j1, j2, j, row_max, row_sum, diag, cur_max)
#endif //BUILD_STRENGTH_MATRIX_PARALLEL
      for (i = 0; i < n; ++i)
        {
          j1 = rows[i];
          j2 = rows[i + 1];
          row_max = 0;
          row_sum = 0;
          diag = data[j1];//data[diag_ind[i]];
          for (j = j1 + 1; j < j2; ++j)
            {
              t_double data_j = data[j];
              t_double row_max_[] = {row_max, data_j};
              t_double row_sum_[] = {0, data_j};
              row_max = row_max_ [data_j < row_max];
              row_sum += row_sum_ [data_j < 0];
            }
#ifdef DEBUG
          //BOSOUT (section::amg, level::debug) << "Row max for row " << i << ":" << row_max << bs_end;
#endif //DEBUG

          if ( -row_sum < diag * max_row_sum)
            continue;

          cur_max = 0;
          for (j = j1 + 1; j < j2; ++j)
            {
              if (data[j] < strength_threshold * row_max)
                {
                  s_markers[j] = 1;
                  cur_max++;
                }
#ifdef DEBUG
              //BOSOUT (section::amg, level::debug) << "S[" << i << ", " << cols[j] << "] = " << s_markers[j] << bs_end;
#endif //DEBUG
            }
#ifdef BUILD_STRENGTH_MATRIX_PARALLEL
          if (cur_max > thread_max)
            {
              thread_max = cur_max;
            }
#else
          if (cur_max > max_connections)
            {
              max_connections = cur_max;
            }
#endif //BUILD_STRENGTH_MATRIX_PARALLEL
        }

#ifdef BUILD_STRENGTH_MATRIX_PARALLEL
      max_connections_local[0] = thread_max;
    }
    for (i = 0; i < n_threads; i++)
      {
        if (max_connections_local[i] > max_connections)
          max_connections = max_connections_local[i];
      }
#endif //BUILD_STRENGTH_MATRIX_PARALLEL



#ifdef DEBUG
    //BOSOUT (section::amg, level::debug) << "*** OUT: build_strength_matrix" << bs_end;
#endif //DEBUG
    return 0;
  }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (simple_smbuilder);
  BLUE_SKY_TYPE_STD_COPY (simple_smbuilder);

  BLUE_SKY_TYPE_IMPL (simple_smbuilder, amg_smbuilder_iface, "simple_smbuilder", "Simple strength matrix builder class", "Realization of simple strength matrix builder");
}  // blue_sky namespace
