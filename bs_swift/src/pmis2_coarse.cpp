/** 
 * @file pmis2_coarse.cpp
 * @brief 
 * @author 
 * @version 
 * @date 2010-03-09
 */

#include "pmis2_coarse.h"
#include "coarse_tools.h"


namespace blue_sky
{
  pmis2_coarse::pmis2_coarse (bs_type_ctor_param) 
                : amg_coarse_iface (),
                  sp_graph (BS_KERNEL.create_object (v_long::bs_type ()))
    {
    }
  pmis2_coarse::pmis2_coarse (const this_t & /*src*/) : bs_refcounter (),
                sp_graph (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }
  
  t_long pmis2_coarse::build (sp_bcsr_t matrix, 
                              spv_double sp_meassure_array, 
                              spv_long sp_cf_markers,
                              spv_long sp_s_markers)
    {
    t_long n = 0;
    t_long n_non_zeros = 0;
    t_long i, j1, j2, j, cl, k;
    t_long f_counter = 0, c_counter = 0;
#ifdef COARSE_PMIS_2_PARALLEL
    int thread_num, n_threads;
    n_threads = omp_get_max_threads ();

#endif //COARSE_PMIS_2_PARALLEL
    t_long last_in_graph = 0;
    t_long ff = -2;
    t_double d;

    // internal checks

    n = matrix->get_n_rows ();
    n_non_zeros = matrix->get_n_non_zeros ();

    t_long *rows = &(*(matrix->get_rows_ptr ()))[0];
    t_long *cols = &(*(matrix->get_cols_ind ()))[0];

    sp_meassure_array->resize (n);
    sp_meassure_array->assign (0);
    sp_graph->resize (n);
    sp_graph->assign (0);

    t_long  *graph = &(*sp_graph)[0];
    t_double *meassure_array = &(*sp_meassure_array)[0];
    t_long  *cf_markers = &(*sp_cf_markers)[0];
    t_long  *s_markers = &(*sp_s_markers)[0];
#ifdef COARSE_PMIS_2_PARALLEL
    assign (last_in_graph_thread, n_threads,0);
#endif //COARSE_PMIS_2_PARALLEL

    // fill meassure array
    for (i = 0; i < n_non_zeros; ++i)
      {
        if (s_markers[i])
          {
            BS_ASSERT (cols[i] < n) (cols[i]) (n);
            meassure_array[cols[i]] += 1;
          }
      }

    // TODO: BUG
    srand (1);

    // TODO: BUG
    // add random number (0,1) to avery node of meassure
    d = (t_double)(1.0 / RAND_MAX);
    for (i = 0, last_in_graph = 0; i < n; ++i)
      {
        if (meassure_array[i] < 1)
          // initial F points not depend from the other points
          {
            cf_markers[i] = -1;
            ++f_counter;
//for comparing results with parallel version (to get equal set of random numbers)
#ifdef SEQ_RAND
            rand();
#endif
          }
        else
          {
            graph[last_in_graph] = i;
            ++last_in_graph;
            cf_markers[i] = ff;
            meassure_array[i] += (t_double)rand () * d;
          }
      }

#ifdef COARSE_PMIS_2_PARALLEL
    t_long *p_last_in_graph_thread=&last_in_graph_thread[0];
#pragma omp parallel private (i, k, j, j1, j2, cl) \
                     firstprivate (ff, graph, last_in_graph_thread) reduction (+: c_counter)

    {
      thread_num = omp_get_thread_num();
      graph += (thread_num * last_in_graph) / n_threads;
      p_last_in_graph_thread += thread_num;
      // IMPORTANT!!! DO NOT SIMPLIFY STATEMENT BELOW!!!
      // It`s not so simple as you might think...
      *p_last_in_graph_thread = ((thread_num + 1) * last_in_graph) / n_threads - (thread_num * last_in_graph) / n_threads;
      //printf ("I am %d graph[0]=%d, last_in_graph_thread = %d of %d\n", thread_num, graph[0], *last_in_graph_thread, last_in_graph);
#else
      t_long *p_last_in_graph_thread = &last_in_graph;
#endif //COARSE_PMIS_2_PARALLEL

      // main loop
      while (last_in_graph > 0)
        {
          // select C points
          coarse_tools::build_independent_set (matrix, 
                                                sp_graph, 
                                                *p_last_in_graph_thread,
                                                sp_meassure_array,
                                                ff,
                                                sp_s_markers,
                                                sp_cf_markers);
          //build_independent_set (matrix, graph, *p_last_in_graph_thread, meassure_array, ff, cf_markers, s_markers);
#ifdef COARSE_PMIS_2_PARALLEL
#pragma omp barrier
#endif //COARSE_PMIS_2_PARALLEL
          // select F points and init C points
          for (k = 0; k < *p_last_in_graph_thread; ++k)
            {
              i = graph[k];
              if (cf_markers[i] == ff)
                // i is new C point
                {
                  cf_markers[i] = (++c_counter);
                  // update graph
                  graph[k--] = graph[--*p_last_in_graph_thread];
                }
              else
                // i is potential new F point
                {
                  j1 = rows[i];
                  j2 = rows[i + 1];
                  for (j = j1; j < j2; ++j)
                    {
                      if (s_markers[j])
                        {
                          cl = cols[j];
                          if (cf_markers[cl] > 0 || cf_markers[cl] == ff)
                            // i is new F point
                            {
                              cf_markers[i] = -1;
                              ++f_counter;
                              graph[k--] = graph[--*p_last_in_graph_thread];
                              break;
                            }
                        }
                    }
                }
            }
          --ff;
#ifdef COARSE_PMIS_2_PARALLEL
#pragma omp master
          {
            last_in_graph = 0;
            for (i = 0; i < n_threads; i++)
              last_in_graph += p_last_in_graph_thread[i];
          }
#pragma omp barrier
        }
#endif //COARSE_PMIS_2_PARALLE

    }
    return c_counter;
    }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (pmis2_coarse);
  BLUE_SKY_TYPE_STD_COPY (pmis2_coarse);

  BLUE_SKY_TYPE_IMPL (pmis2_coarse, amg_coarse_iface, "pmis2_coarse", "PMIS coarsering builder class", "Realization of PMIS coarsering");
}  // blue_sky namespace
