/** 
 * @file cljp_coarse.cpp
 * @brief implementation of Cleary-Luby-Jones-Plassmann grid coasering method
 * @author 
 * @version 
 * @date 2010-03-10
 */
#include "cljp_coarse.h"

#include "coarse_tools.h"


namespace blue_sky
{
  cljp_coarse::cljp_coarse (bs_type_ctor_param) 
               : amg_coarse_iface (),
                 sp_graph (BS_KERNEL.create_object (v_long::bs_type ()))
    {
    }
  cljp_coarse::cljp_coarse (const this_t & /*src*/) : bs_refcounter (),
               sp_graph (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }
  
  t_long 
  cljp_coarse::build (sp_bcsr_t matrix, 
                      spv_double sp_meassure_array, 
                      spv_long sp_cf_markers,
                      spv_long sp_s_markers)
    {
      t_long n = 0;
      t_long n_non_zeros = 0;
      t_long i, j1, j2, j, jj, jj1, jj2, cl, cl2, k;
      t_double d;
      t_long last_in_graph = 0;
      t_long ff = -2;
      t_long iter;
      const t_long common_c_pnt = 2;
      t_long c_counter = 0;

      n = matrix->get_n_rows ();
      n_non_zeros = matrix->get_n_non_zeros ();

      t_long *rows = &(*(matrix->get_rows_ptr ()))[0];
      t_long *cols = &(*(matrix->get_cols_ind ()))[0];

      sp_graph->resize (n);
      sp_graph->assign (0);
      sp_meassure_array->resize (n);
      sp_meassure_array->assign (0);
      sp_cf_markers->resize (n);
      sp_cf_markers->assign (0);

      t_long  *graph = &(*sp_graph)[0];
      t_double *meassure_array = &(*sp_meassure_array)[0];
      t_long  *cf_markers = &(*sp_cf_markers)[0];
      t_long  *s_markers = &(*sp_s_markers)[0];
      

      // fill meassure array
      for (i = 0; i < n_non_zeros; ++i)
        {
          if (s_markers[i])
            {
              BS_ASSERT (cols[i] < n);
              meassure_array[cols[i]] += 1.0;
            }
        }

      // add random number to meassure array
      d = (t_double)(1.0 / RAND_MAX);
      //TODO: Ask Oleg! about srand
      srand (66);
#ifdef AMG_COARSE_CLJP_PARALLEL
#pragma omp parallel for
#endif //AMG_COARSE_CLJP_PARALLEL
      for (i = 0; i < n; ++i)
        meassure_array[i] += (t_double)rand () * d;

      for (i = 0, last_in_graph = 0; i < n; ++i)
        {
          graph[last_in_graph++] = i;
          cf_markers[i] = ff;
        }

      // main loop
      for (iter = 0; last_in_graph > 0; ++iter)
        {
          coarse_tools::build_independent_set (matrix, 
                                               sp_graph, 
                                               last_in_graph,
                                               sp_meassure_array,
                                               ff,
                                               sp_s_markers,
                                               sp_cf_markers);
          //build_independent_set (matrix, graph, last_in_graph, meassure_array, ff, cf_markers, s_markers);

          // update weights (meassure_array)
#ifdef AMG_COARSE_CLJP_PARALLEL
#pragma omp parallel for private (i, j1, j2, j, cl, jj1, jj2, jj, cl2) \
                             reduction (+ : c_counter)
#endif //AMG_COARSE_CLJP_PARALLEL
          for (k = 0; k < last_in_graph; ++k)
            {
              i = graph[k];
              if (cf_markers[i] == ff || cf_markers[i] > 0)
                {
                  cf_markers[i] = 1;
                  j1 = rows[i];
                  j2 = rows[i + 1];
                  for (j = j1; j < j2; ++j)
                    {
                      if (s_markers[j] > 0)
                        {
#ifdef AMG_COARSE_CLJP_PARALLEL
#pragma omp critical (cljp)
                          {
#endif //AMG_COARSE_CLJP_PARALLEL
                            s_markers[j] = -1;
                            --meassure_array[cols[j]];
#ifdef AMG_COARSE_CLJP_PARALLEL
                          }
#endif //AMG_COARSE_CLJP_PARALLEL
                        }
                    }
                  ++c_counter;
                }
              else
                {
                  j1 = rows[i];
                  j2 = rows[i + 1];
                  for (j = j1; j < j2; ++j)
                    {
                      if (s_markers[j])
                        {
                          cl = cols[j];
                          if (cf_markers[cl] == ff || cf_markers[cl] > 0)
                            {
#ifdef AMG_COARSE_CLJP_PARALLEL
#pragma omp critical (cljp)
                              {
#endif //AMG_COARSE_CLJP_PARALLEL
                                s_markers[j] = -1;
                                cf_markers[cl] = common_c_pnt;
#ifdef AMG_COARSE_CLJP_PARALLEL
                              }
#endif //AMG_COARSE_CLJP_PARALLEL
                            }
                        }
                    }
                  for (j = j1; j < j2; ++j)
                    {
                      if (s_markers[j] > 0)
                        {
                          cl = cols[j];
                          // check for common C point
                          jj1 = rows[cl];
                          jj2 = rows[cl + 1];
                          for (jj = jj1; jj < jj2; ++jj)
                            {
                              if (s_markers[jj])
                                {
                                  cl2 = cols[jj];
                                  if (cf_markers[cl2] == common_c_pnt)
                                    {
#ifdef AMG_COARSE_CLJP_PARALLEL
#pragma omp critical (cljp)
                                      {
#endif //AMG_COARSE_CLJP_PARALLEL
                                        s_markers[j] = -1;
                                        --meassure_array[cl];
#ifdef AMG_COARSE_CLJP_PARALLEL
                                      }
#endif //AMG_COARSE_CLJP_PARALLEL
                                      break;
                                    }
                                }
                            }
                        }
                    }
                  for (j = j1; j < j2; ++j)
                    {
                      cl = cols[j];
                      if (cf_markers[cl] == common_c_pnt)
#ifdef AMG_COARSE_CLJP_PARALLEL
#pragma omp critical (cljp)
#endif //AMG_COARSE_CLJP_PARALLEL
                        cf_markers[cl] = 1;

                    }
                }
            }
          // find new F points and update graph

          for (k = 0; k < last_in_graph; ++k)
            {
              i = graph[k];
              if (cf_markers[i] > 0 || meassure_array[i] < 1)
                {
                  graph[k] = graph[last_in_graph - 1];
                  --last_in_graph;
                  --k;
                }
            }
          --ff;
        }
      return c_counter;
    }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (cljp_coarse);
  BLUE_SKY_TYPE_STD_COPY (cljp_coarse);

  BLUE_SKY_TYPE_IMPL (cljp_coarse, amg_coarse_iface, "cljp_coarse", "CLJP coarsering builder class", "Realization of CLJP coarsering");
}  // blue_sky namespace
