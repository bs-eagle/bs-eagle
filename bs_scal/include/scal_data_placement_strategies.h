/**
 * \file scal_data_placement_strategies.h
 * \brief different strategies of scal data placement
 * \author Sergey Miryanov
 * \date 22.05.2008
 * */
#ifndef BS_SCAL_DATA_PLACEMENT_STRATEGIES_H_
#define BS_SCAL_DATA_PLACEMENT_STRATEGIES_H_

#include "scal_3p.h"
#include "scal_data_source.h"
#include "scal_data_vector.h"
#include "scal_data_placement_info.h"

namespace blue_sky
  {
  namespace scal
    {
    namespace data_placement
      {

      struct separete_vectors_t
        {
          static double *
          place_data (double *memory, scal::data_placement::scal_placement_info &info, const input_scal_region &region)
          {
            // copy data
            memory = copy (memory, region.Sp);
            memory = copy (memory, region.So);
            memory = copy (memory, region.Krp);
            memory = copy (memory, region.Krop);
            memory = copy (memory, region.Pcp);

            // set info
            info.sp_step      = 1;
            info.so_step      = 1;
            info.krp_step     = 1;
            info.krop_step    = 1;
            info.pcp_step     = 1;
            info.sp_offset    = 0;
            info.so_offset    = info.sp_offset + (int)region.Sp.size ();
            info.krp_offset   = info.so_offset + (int)region.So.size ();
            info.krop_offset  = info.krp_offset + (int)region.Krp.size ();
            info.pcp_offset   = info.krop_offset + (int)region.Krop.size ();

            return memory;
          }

private:

          static double *
          copy (double *memory, const input_scal_region::vector_t &src)
          {
            size_t size = sizeof (input_scal_region::vector_t::value_type) * src.size ();
            memcpy (memory, &src[0], size);

            return memory + src.size ();
          }

        };

      struct spkrp_sokrop_pcp_t
        {
          static double *
          place_data (double *memory, scal::data_placement::scal_placement_info &info, const input_scal_region &region)
          {
            // copy data
            for (int i = 0, cnt = (int)region.Sp.size (); i < cnt; ++i, memory += 2)
              {
                memory[0] = region.Sp[i];
                memory[1] = region.Krp[i];
              }
            for (int i = 0, cnt = (int)region.So.size (); i < cnt; ++i, memory += 2)
              {
                memory[0] = region.So[i];
                memory[1] = region.Krop[i];
              }
            for (int i = 0, cnt = (int)region.Pcp.size (); i < cnt; ++i, memory += 1)
              memory[0] = region.Pcp[i];

            // set info
            info.sp_step      = 2;
            info.so_step      = 2;
            info.krp_step     = 2;
            info.krop_step    = 2;
            info.pcp_step     = 1;
            info.sp_offset    = 0;
            info.krp_offset   = 1;
            info.so_offset    = info.sp_step    + info.sp_step * (int)region.Sp.size ();
            info.krop_offset  = info.krp_offset + info.sp_step * (int)region.Sp.size ();
            info.pcp_offset   = info.so_offset  + info.so_step * (int)region.So.size ();

            return memory;
          }
        };

      struct all_regions_t
        {
          typedef scal_2p_data_holder                           scal_2p_data_holder_t;
          typedef v_double                                      array_item_t;
          typedef smart_ptr <v_double, true>                    sp_array_item_t;
          
          static void
          place_spof_data (sp_array_item_t dst, scal::data_placement::scal_placement_info &info, const sp_array_item_t src, bool is_water)
          {
            t_double *dst_array       = dst->data ();
            t_double const *src_array = src->data ();
            size_t size = dst->size () - (5 * (src->size () / 4));
            size_t src_size = src->size ();
            for (size_t i = 0, j = 0, cnt = src->size (); i < cnt; i += 4, j += 5)
              {
                dst_array[size + j + 0] = src_array[i + 0];
                dst_array[size + j + 1] = 1.0 - src_array[src_size - (i + 4) + 0];
                dst_array[size + j + 2] = src_array[i + 1];
                dst_array[size + j + 3] = src_array[src_size - (i + 4) + 2];
                dst_array[size + j + 4] = is_water ? -src_array[i + 3] : src_array[i + 3];
              }

            info.sp_step      = 5;
            info.so_step      = 5;
            info.krp_step     = 5;
            info.krop_step    = 5;
            info.pcp_step     = 5;
            info.sp_offset    = 0;
            info.so_offset    = 1;
            info.krp_offset   = 2;
            info.krop_offset  = 3;
            info.pcp_offset   = 4;
          }

          static void
          place_spfn_data (sp_array_item_t dst, scal_placement_info &info, const sp_array_item_t src, bool is_water)
          {
            array_item_t &dst_array = *dst;
            array_item_t &src_array = *src;
            size_t size = dst->size () - src->size ();
            for (size_t i = 0, cnt = src->size (); i < cnt; i += 3)
              {
                dst_array[size + i + 0] = src_array[i + 0];
                dst_array[size + i + 1] = src_array[i + 1];
                dst_array[size + i + 2] = is_water ? -src_array[i + 2] : src_array[i + 2];
              }

            info.sp_step      = 3;
            info.krp_step     = 3;
            info.pcp_step     = 3;
            info.so_step      = 2;
            info.krop_step    = 2;
            info.sp_offset    = 0;
            info.krp_offset   = 1;
            info.pcp_offset   = 2;
            info.so_offset    = 0;
            info.krop_offset  = 1;
          }

          static void
          place_sof3_data (sp_array_item_t dst, scal_placement_info &info, const sp_array_item_t src, bool is_water)
          {
            array_item_t &dst_array = *dst;
            array_item_t &src_array = *src;
            size_t size = dst->size () - src->size ();
            for (size_t i = 0, j = 0, cnt = src->size (); i < cnt; i += 3, j += 2)
              {
                dst_array[size + j + 0] = src_array[i + 0];
                dst_array[size + j + 1] = src_array[i + (is_water ? 1 : 2)];
              }

            info.sp_step      = 3;
            info.krp_step     = 3;
            info.pcp_step     = 3;
            info.so_step      = 2;
            info.krop_step    = 2;
            info.sp_offset    = 0;
            info.krp_offset   = 1;
            info.pcp_offset   = 2;
            info.so_offset    = 0;
            info.krop_offset  = 1;
          }
        };

    } // namespace data_placement
  } // namespace scal
} // namespace blue_sky


#endif // #ifndef BS_SCAL_DATA_PLACEMENT_STRATEGIES_H_
