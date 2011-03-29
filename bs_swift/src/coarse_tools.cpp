/** 
 * @file coarse_tools.cpp
 * @brief 
 * @author 
 * @version 
 * @date 2010-03-03
 */

#include "coarse_tools.h"

namespace blue_sky
  {

  coarse_tools::coarse_tools (bs_type_ctor_param) 
    {
    }
  coarse_tools::coarse_tools (const this_t & /*src*/) : bs_refcounter (), objbase ()
     {
     }

//build transposed strength matrix
int 
coarse_tools::build_transposed_strength_matrix (sp_bcsr_t matrix, 
                                                spv_long sp_s_markers,
                                                spv_long sp_st_markers,
                                                spv_long sp_st_cols,
                                                spv_long sp_st_rows)
  {
    t_long cl, i, j1, j2, j;

    t_long nnz = matrix->get_n_non_zeros ();
    t_long n = matrix->get_n_rows ();
    if (nnz < 1 || n < 1)
      return -88;

    t_long *rows = &(*(matrix->get_rows_ptr ()))[0];
    t_long *cols = &(*(matrix->get_cols_ind ()))[0];

    sp_st_markers->resize (nnz);
    sp_st_markers->assign (0);

    sp_st_cols->resize (nnz);
    sp_st_cols->assign (0);

    sp_st_rows->resize (n + 1);
    t_long  *st_rows = &(*sp_st_rows)[0];
    memcpy (st_rows, rows, (n + 1) * sizeof (t_long));

    t_long  *s_markers = &(*sp_s_markers)[0];
    t_long  *st_markers = &(*sp_st_markers)[0];
    t_long  *st_cols = &(*sp_st_cols)[0];


    // build S^T matrix
    for (i = 0; i < n; ++i)
      {
        j1 = rows[i];
        j2 = rows[i + 1];
        for (j = j1; j < j2; ++j)
          {
            cl = cols[j];
            BS_ASSERT (st_rows[cl] < nnz) (st_rows[cl]) (nnz);
            BS_ASSERT (st_rows[cl] < nnz) (st_rows[cl]) (nnz);

            st_markers[st_rows[cl]] += s_markers[j];
            st_cols[st_rows[cl]] = i;
            ++st_rows[cl];
          }
      }

    return 0;
  }

 int 
coarse_tools::build_independent_set (sp_bcsr_t matrix,
                                        spv_long sp_graph,
                                        const t_long last_in_graph,
                                        spv_double sp_meassure_array,
                                        const t_long ff,
                                        spv_long sp_s_markers,
                                        spv_long sp_cf_markers)
  {
    t_long k, i, j, j1, j2, cl;
    t_double cur_m;

    t_long *rows = &(*(matrix->get_rows_ptr ()))[0];
    t_long *cols = &(*(matrix->get_cols_ind ()))[0];
    t_long  *graph = &(*sp_graph)[0];
    t_double *meassure_array = &(*sp_meassure_array)[0];
    t_long  *cf_markers = &(*sp_cf_markers)[0];
    t_long  *s_markers = &(*sp_s_markers)[0];

#ifdef AMG_COARSE_CLJP_PARALLEL
#pragma omp parallel for private (i, cur_m, j1, j2, j, cl)
#endif //AMG_COARSE_CLJP_PARALLEL
    for (k = 0; k < last_in_graph; ++k)
      {
        i = graph[k];
        cur_m = meassure_array[i];
        if (cur_m > 1)
          {
            j1 = rows[i];
            j2 = rows[i + 1];

            for (j = j1; j < j2; ++j)
              {
                if (s_markers[j])
                  {
                    cl = cols[j];
                    // if i strongly depends from cl
                    if (cur_m > meassure_array[cl] && cf_markers[cl] == ff)
                      // cl not a C point move it to the next iteration
                      cf_markers[cl] = ff - 1;
                    //else if (cur_m < meassure_array[cl] && cf_markers[cl] == ff)
                    //  cf_markers[i] = ff - 1;
                  }
              }
          }
      }
    return 0;
  }
  BLUE_SKY_TYPE_STD_CREATE (coarse_tools);
  BLUE_SKY_TYPE_STD_COPY (coarse_tools);

  BLUE_SKY_TYPE_IMPL (coarse_tools, objbase, "coarse_tools", "Coarsering tools", "Realization of coarsering tools");
  };
