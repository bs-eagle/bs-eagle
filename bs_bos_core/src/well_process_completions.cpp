/**
 *       \file  wells_process_completions.cpp
 *      \brief  process connections for well
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.07.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "well_process_completions.h"
#include "calc_well.h"
#include "wells_common_const.h"
#include "fi_params.h"
#include "well_connection.h"
#include "well_pool_iface.h"
//#include "frac_comp_ident.h"

namespace blue_sky
  {
  namespace wells
    {
      void
      completion::process_completion (const physical_constants &internal_constants,
                                      const sp_params_t &params,
                                      const sp_mesh_iface_t &mesh,
                                      const stdv_float &perm,
                                      const stdv_float &ntg,
                                      sp_well_t &well,
                                      std::vector <well_hit_cell_3d> &well_path_segs,
                                      t_uint /*con_branch*/,
                                      t_float con_md,
                                      t_float con_len)
        {
          t_float well_seg_start_point[3];
          t_float well_seg_end_point[3];
          const smart_ptr <rs_smesh_iface, true> s_mesh(mesh, bs_dynamic_cast ());

          // loop through well segments to find segment containing completion
          for (t_ulong i = 0; i < well_path_segs.size () - 1; ++i)
            {
              ulong cell_num  = well_path_segs[i].cell;
              t_float well_md = well_path_segs[i].md;
              well_seg_start_point[0] = well_path_segs[i].where.x ();
              well_seg_start_point[1] = well_path_segs[i].where.y ();
              well_seg_start_point[2] = well_path_segs[i].where.z ();
              //uint facet      = well_path_segs[i].facet;
              //bool is_node    = well_path_segs[i].is_node;

              t_float well_md_next = 0;
              t_float well_length = 0;

              well_seg_end_point[0] = well_path_segs[i + 1].where.x ();
              well_seg_end_point[1] = well_path_segs[i + 1].where.y ();
              well_seg_end_point[2] = well_path_segs[i + 1].where.z ();
              well_md_next = well_path_segs[i + 1].md;
              well_length = well_md_next - well_md;

              // check completion position with well path
              if ((con_md > well_md + well_length) || (con_md + con_len < well_md))
                continue;

              t_long i_coord, j_coord, k_coord;
              completion_coords_t completion_data;
              s_mesh->get_element_int_to_ijk (cell_num, i_coord, j_coord, k_coord);

              // completion fully inside well segment
              if (con_md > well_md && (con_md + con_len) < (well_md + well_length))
                {
                  for (t_uint j = 0; j < 3; ++j)
                    {
                      completion_data.x1[j] = well_seg_start_point[j] + (con_md - well_md) / well_length * (well_seg_end_point[j] - well_seg_start_point[j]);
                      completion_data.x2[j] = well_seg_start_point[j] + (con_md + con_len - well_md) / well_length * (well_seg_end_point[j] - well_seg_start_point[j]);
                    }

                }
              // completion occupy all well segment
              else if (con_md <= well_md && (con_md + con_len) >= (well_md + well_length))
                {
                  for (t_uint j = 0; j < 3; ++j)
                    {
                      completion_data.x1[j] = well_seg_start_point[j];
                      completion_data.x2[j] = well_seg_end_point[j];
                    }
                }
              else if (con_md >= well_md && (con_md + con_len) >= (well_md + well_length))
                {
                  for (t_uint j = 0; j < 3; ++j)
                    {
                      completion_data.x1[j] = well_seg_start_point[j] + (con_md - well_md) / well_length * (well_seg_end_point[j] - well_seg_start_point[j]);
                      completion_data.x2[j] = well_seg_end_point[j];
                    }
                }
              else if (con_md <= well_md && (con_md + con_len) <= (well_md + well_length))
                {
                  for (t_uint j = 0; j < 3; ++j)
                    {
                      completion_data.x1[j] = well_seg_start_point[j];
                      completion_data.x2[j] = well_seg_start_point[j] + (con_md + con_len - well_md) / well_length * (well_seg_end_point[j] - well_seg_start_point[j]);
                    }
                }
              else
                {
                  bs_throw_exception ("Error: undefined completion position along well path");
                }

              const sp_connection_t &con = well->get_connection_map (cell_num);
              if (!con)
                {
                  const sp_connection_t con1 = well->add_completion (i_coord, j_coord, k_coord, cell_num, completion_data);
                  //con->compute_factors (internal_constants, params, mesh, perm, ntg, true, true);
                  compute_factors::completion_connection_factor::compute (*con1, internal_constants, params, mesh, perm, ntg);
                }
              else
                {
                  BS_ASSERT (con->get_CFF_flag () == true);
                  well_t::sp_connection_t temp_con = BS_KERNEL.create_object (well_t::connection_t::bs_type (), true);
                  temp_con->set_coord (i_coord, j_coord, k_coord, cell_num);
                  temp_con->set_completion_coords (completion_data);
                  //temp_con->compute_factors (internal_constants, params, mesh, perm, ntg, true, true);
                  compute_factors::completion_connection_factor::compute (*temp_con, internal_constants, params, mesh, perm, ntg);
                  // set new factor as sum of all connection factors
                  con->set_factor (con->get_fact () + temp_con->get_fact ());
                }

            }

        }



  } // namespace wells
} // namespace blue_sky
