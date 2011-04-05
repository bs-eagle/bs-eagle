/** 
 * @file ruge_coarse.cpp
 * @brief implimentation of Ruge algorithm
 * @author 
 * @version 
 * @date 2010-03-03
 */
#include "ruge_coarse.h"
#include "coarse_tools.h"
#include "link_list.h"


namespace blue_sky
{
  ruge_coarse::ruge_coarse (bs_type_ctor_param) 
                : amg_coarse_iface (),
                sp_st_rows (BS_KERNEL.create_object (v_long::bs_type ())),
                sp_st_cols (BS_KERNEL.create_object (v_long::bs_type ())),
                sp_st_markers (BS_KERNEL.create_object (v_long::bs_type ())),
                sp_next (BS_KERNEL.create_object (v_long::bs_type ())),
                sp_prev (BS_KERNEL.create_object (v_long::bs_type ()))
    {
    }
  
  ruge_coarse::ruge_coarse (const this_t & /*src*/) : bs_refcounter (),
                        sp_st_rows (BS_KERNEL.create_object (v_long::bs_type ())),
                        sp_st_cols (BS_KERNEL.create_object (v_long::bs_type ())),
                        sp_st_markers (BS_KERNEL.create_object (v_long::bs_type ())),
                        sp_next (BS_KERNEL.create_object (v_long::bs_type ())),
                        sp_prev (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }
  
#define C_PT 1
#define F_PT -1
#define Z_PT -2
#define SF_PT -3  /* special fine points */
#define UNDECIDED 0
   
  t_long ruge_coarse::build (sp_bcsr_t matrix, 
                              spv_double sp_meassure_array, 
                              spv_long sp_cf_markers,
                              spv_long sp_s_markers)
    {
      t_long n = 0;
      t_long n_non_zeros = 0;
      t_long i, j1, j2, j, jj, jj1, jj2;
      t_long num_left, c_counter = 0;
      t_long meassure, new_meas;
      t_long elem, elem2, index;
      t_long row_start = 0, row_end;

      double_linked_list   *LoL_head;
      double_linked_list   *LoL_tail;

#ifdef AMG_COARSE_RUGE_PARALLEL
      int thread_num, n_threads = 1;
      n_threads = omp_get_max_threads ();

#endif //AMG_COARSE_RUGE_PARALLEL

      // internal checks
      //BS_ASSERT (matrix);
      //BS_ASSERT (s_markers.size ());

      // prepare ST matrix
      //build_transposed_strength_matrix (matrix, n_memory_st_markers,st_markers, st_cols, s_markers);
      coarse_tools::build_transposed_strength_matrix (matrix, sp_s_markers, sp_st_markers, sp_st_cols, sp_st_rows);

      n = matrix->get_n_rows ();
      row_end = n;
      n_non_zeros = matrix->get_n_non_zeros ();

      t_long *rows = &(*(matrix->get_rows_ptr ()))[0];
      t_long *cols = &(*(matrix->get_cols_ind ()))[0];

      sp_next->resize (n);
      sp_next->assign (0);
      sp_prev->resize (n);
      sp_prev->assign (0);

      sp_meassure_array->resize (n);
      sp_meassure_array->assign (0);

      LoL_head = 0;
      LoL_tail = 0;

      sp_cf_markers->resize (n);
      sp_cf_markers->assign (0);

      t_double *meassure_array = &(*sp_meassure_array)[0];
      t_long  *cf_markers = &(*sp_cf_markers)[0];
      t_long  *s_markers = &(*sp_s_markers)[0];
      t_long  *st_markers = &(*sp_st_markers)[0];
      t_long  *st_cols = &(*sp_st_cols)[0];
      t_long  *prev = &(*sp_prev)[0];
      t_long  *next = &(*sp_next)[0];

      //t_long  *st_rows = &(*sp_st_rows)[0];
      // fill meassure array
      for (i = 0; i < n_non_zeros; ++i)
        {
          if (s_markers[i])
            {
              BS_ASSERT (cols[i] < n);
              meassure_array[cols[i]] += 1;
            }
        }
      num_left = n;

#ifdef AMG_COARSE_RUGE_PARALLEL
#pragma omp parallel private (thread_num, row_start, row_end, num_left, i, meassure, \
                                j, j1, j2, elem, new_meas, index, jj, jj1, jj2, elem2) \
                       firstprivate (LoL_head, LoL_tail) \
                       reduction (+ : c_counter)
      {
        thread_num = omp_get_thread_num ();
        row_end = n * (thread_num + 1) / n_threads;
        row_start = n * thread_num / n_threads;
        num_left = row_end - row_start;
#else //AMG_COARSE_RUGE_PARALLEL
      {
#endif //AMG_COARSE_RUGE_PARALLEL

        for (i = row_start; i < row_end; i++)
          {
            meassure = meassure_array[i];
            if (meassure > 0)
              {
                enter_on_lists (&LoL_head, &LoL_tail, meassure, i, &next[0], &prev[0]);
              }
            else
              {
                cf_markers[i] = F_PT;
                j1 = rows[i];
                j2 = rows[i + 1];
                for (j = j1; j < j2; j++)
                  {
                    if (s_markers[j] > 0)
                      {
                        elem = cols[j];
                        if (elem >= row_start && elem < row_end)
                          {
                            if (elem < i)
                              {
                                t_double new_meas = meassure_array[elem];
                                if (new_meas > 0)
                                  remove_point(&LoL_head, &LoL_tail, new_meas, elem, &next[0], &prev[0]);
                                new_meas = ++(meassure_array[elem]);
                                enter_on_lists(&LoL_head, &LoL_tail, new_meas, elem, &next[0], &prev[0]);
                              }
                            else
                              {
                                ++(meassure_array[elem]);
                              }
                          }
                      }
                  }
                --num_left;
              }

          }

        // main loop
        while (num_left > 0)
          {
            index = LoL_head -> head;

            cf_markers[index] = 1;
            ++c_counter;
            meassure = meassure_array[index];
            meassure_array[index] = 0;
            --num_left;

            remove_point(&LoL_head, &LoL_tail, meassure, index, &next[0], &prev[0]);
            j1 = rows[index];
            j2 = rows[index + 1];

            for (j = j1; j < j2; j++)
              {
                if (st_markers[j] > 0)
                  {
                    elem = st_cols[j];
                    if (cf_markers[elem] == UNDECIDED && elem >= row_start && elem < row_end)
                      {
                        cf_markers[elem] = F_PT;
                        meassure = meassure_array[elem];

                        remove_point(&LoL_head, &LoL_tail, meassure, elem, &next[0], &prev[0]);
                        --num_left;
                        jj1 = rows[elem];
                        jj2 = rows[elem + 1];
                        for (jj = jj1; jj < jj2; jj++)
                          {
                            if (s_markers[jj] > 0)
                              {
                                elem2 = cols[jj];
                                if (cf_markers[elem2] == UNDECIDED && elem2 >= row_start && elem2 < row_end)
                                  {

                                    meassure = meassure_array[elem2];
                                    remove_point(&LoL_head, &LoL_tail, meassure,
                                                 elem2, &next[0], &prev[0]);

                                    new_meas = ++(meassure_array[elem2]);

                                    enter_on_lists(&LoL_head, &LoL_tail, new_meas,
                                                   elem2, &next[0], &prev[0]);
                                  }
                              }
                          }
                      }
                  }
              }

            for (j = j1; j < j2; j++)
              {
                if (s_markers[j] > 0)
                  {
                    elem = cols[j];
                    if (cf_markers[elem] == UNDECIDED && elem >= row_start && elem < row_end)
                      {
                        meassure = meassure_array[elem];

                        remove_point(&LoL_head, &LoL_tail, meassure, elem,  &next[0], &prev[0]);

                        meassure_array[elem] = --meassure;

                        if (meassure > 0)
                          {
                            enter_on_lists(&LoL_head, &LoL_tail, meassure, elem,  &next[0], &prev[0]);
                          }
                        else
                          {
                            cf_markers[elem] = F_PT;
                            --num_left;
                            jj1 = rows[elem];
                            jj2 = rows[elem + 1];
                            for (jj = jj1; jj < jj2; jj++)
                              {
                                if (s_markers[jj] > 0)
                                  {
                                    elem2 = cols[jj];
                                    if (cf_markers[elem2] == UNDECIDED && elem2 >= row_start && elem2 < row_end)
                                      {
                                        new_meas = meassure_array[elem2];
                                        remove_point(&LoL_head, &LoL_tail, new_meas, elem2,  &next[0], &prev[0]);
                                        new_meas = ++(meassure_array[elem2]);
                                        enter_on_lists(&LoL_head, &LoL_tail, new_meas, elem2,  &next[0], &prev[0]);
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
          }
      } //end parallel
      return c_counter;
    }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (ruge_coarse);
  BLUE_SKY_TYPE_STD_COPY (ruge_coarse);

  BLUE_SKY_TYPE_IMPL (ruge_coarse, amg_coarse_iface, "ruge_coarse", "Ruge coarsering builder class", "Realization of Ruge coarsering");
}  // blue_sky namespace

