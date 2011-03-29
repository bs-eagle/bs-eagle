/** 
 * @file pmis_coarse.cpp
 * @brief implimentation of PMIS algorithm
 * @author 
 * @version 
 * @date 2010-03-03
 */
#include "pmis_coarse.h"
#include "coarse_tools.h"


namespace blue_sky
{
  pmis_coarse::pmis_coarse (bs_type_ctor_param) 
                : amg_coarse_iface (),
                  sp_st_rows (BS_KERNEL.create_object (v_long::bs_type ())),
                  sp_st_cols (BS_KERNEL.create_object (v_long::bs_type ())),
                  sp_st_markers (BS_KERNEL.create_object (v_long::bs_type ()))
    {
    }
  pmis_coarse::pmis_coarse (const this_t & /*src*/) : bs_refcounter (),
                sp_st_rows (BS_KERNEL.create_object (v_long::bs_type ())),
                sp_st_cols (BS_KERNEL.create_object (v_long::bs_type ())),
                sp_st_markers (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }
  
  t_long pmis_coarse::build (sp_bcsr_t a_matrix, 
                                      spv_double sp_meassure_array, 
                                      spv_long sp_cf_markers,
                                      spv_long sp_s_markers)
    {
      t_long n = 0;
      t_long n_non_zeros = 0;
      t_long i, j1, j2, j, cl;
      t_double d;
      t_long f_counter = 0, c_counter = 0;
      t_double cur_meassure;
      t_long b_ind = 0;

      n = a_matrix->get_n_rows ();
      n_non_zeros = a_matrix->get_n_non_zeros ();

      t_long *rows = &(*(a_matrix->get_rows_ptr ()))[0];
      t_long *cols = &(*(a_matrix->get_cols_ind ()))[0];

      // initialize cf_markers
      sp_meassure_array->resize (n);
      sp_meassure_array->assign (0);
      sp_cf_markers->resize (n);
      sp_cf_markers->assign (0);

      // prepare ST matrix
      coarse_tools::build_transposed_strength_matrix (a_matrix, sp_s_markers, sp_st_markers, sp_st_cols, sp_st_rows);

      t_double *meassure_array = &(*sp_meassure_array)[0];
      t_long  *cf_markers = &(*sp_cf_markers)[0];
      t_long  *s_markers = &(*sp_s_markers)[0];
      t_long  *st_markers = &(*sp_st_markers)[0];
      t_long  *st_cols = &(*sp_st_cols)[0];
      //t_long  *st_rows = &(*sp_st_rows)[0];

      // fill meassure array
      for (i = 0; i < n_non_zeros; ++i)
        {
          if (s_markers[i])
            {
              BS_ASSERT (cols[i] < n);
              meassure_array[cols[i]] += 1.0;
            }
        }
      // add random number (0,1) to avery node of meassure
      d = (t_double)(1.0 / RAND_MAX);

      for (i = 0; i < n; ++i)
        {
          if (meassure_array[i] < 1)
            // initial F points not depend from the other points
            cf_markers[i] = - (++f_counter);
          else
            meassure_array[i] += (t_double)rand () * d;
        }

      // main loop
      while (f_counter + c_counter < n)
        {
          // find first not market point
          for (i = b_ind; i < n && cf_markers[i]; ++i)
            ;
          b_ind = i;
          if (i >= n)
            return -1;
          cur_meassure = meassure_array[i];
          // find point with maximum weight
          for (;;)
            {
              j1 = rows[i];
              j2 = rows[i + 1];
              for (j = j1; j < j2; ++j)
                {
                  if (st_markers[j] > 0)
                    {
                      cl = st_cols[j];
                      if (!cf_markers[cl] && cur_meassure < meassure_array[cl])
                        {
                          i = cl;
                          cur_meassure = meassure_array[cl];
                          break;
                        }
                    }
                }
              if (j == j2)
                {
                  // i is the new C point
                  cf_markers[i] = (++c_counter);
                  // mark all dependent F points
                  for (j = j1; j < j2; ++j)
                    {
                      if (s_markers[j] > 0)
                        {
                          cl = cols[j];
                          if (!cf_markers[cl])
                            cf_markers[cl] = -(++f_counter);
                        }
                    }
                  break;
                }
            }
        }
      return c_counter;
    }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (pmis_coarse);
  BLUE_SKY_TYPE_STD_COPY (pmis_coarse);

  BLUE_SKY_TYPE_IMPL (pmis_coarse, amg_coarse_iface, "pmis_coarse", "PMIS coarsering builder class", "Realization of PMIS coarsering");
}  // blue_sky namespace

