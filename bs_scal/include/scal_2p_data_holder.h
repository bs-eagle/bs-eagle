/**
 * \file scal_2p_data_holder.h
 * \brief holder for SCAL data
 * \author Sergey Miryanov
 * \date 30.12.2008
 * */
#ifndef BS_SCAL_SCAL_2P_DATA_HOLDER_H_
#define BS_SCAL_SCAL_2P_DATA_HOLDER_H_

#include "scal_data_placement_info.h"

namespace blue_sky {

  /**
   * \brief hold scal region data (So, Sp, Krp, Krop)
   * \detail we hold data as a plain memory region and a list of region_info
   * */
  class BS_API_PLUGIN scal_2p_data_holder : public scal_2p_data_holder_iface
  {
  public:
    typedef t_double                                item_t;
    typedef t_int                                   index_t;

    typedef scal_2p_data_holder                     this_t;

    typedef scal_region                             scal_region_t;
    typedef scal_region_info <item_t>               scal_region_info_t;
    typedef data_vector <item_t>                    data_vector_t;
    typedef v_double                                item_array_t;
    typedef smart_ptr <item_array_t, true>          sp_array_item_t;
    typedef shared_vector <scal_region_info_t>      region_vector_t;
    typedef shared_vector <scal_region_t *>         region_vector_2_t;
    typedef smart_ptr <this_t, true>								sp_scal_data_t;

    void add_spof (sp_array_item_t const &data, bool is_water);
    void add_spfn (sp_array_item_t const &data, t_long region_index, bool is_water);
    void add_sof3 (sp_array_item_t const &data, t_long region_index, bool is_water);

    ~scal_2p_data_holder ()
    {
      for (size_t i = 0, cnt = region_2_.size (); i < cnt; ++i)
        {
          delete region_2_[i];
        }
    }

    t_float
    get_phase_sat_min (t_long region) const
    {
      return get_region (region).get_phase_sat_min ();
    }

    t_float
    get_phase_sat_max (t_long region) const
    {
      return get_region (region).get_phase_sat_max ();
    }

    t_float
    get_pcp_max (t_long region) const
    {
      return get_region (region).get_pcp_max ();
    }

    const scal_region_t &
    get_region (index_t index) const
    {
      return *region_2_[index];
    }
    scal_region_t &
    get_region (index_t index)
    {
      return *region_2_[index];
    }

    // returns region from region_ not from region_2_ (region_2_ precalculated)
    scal_region_t 
    get_region_from_info (int index) const
    {
      item_array_t &data_array = *data_;
      BS_ASSERT (index >= 0 && index < (int)region_.size ());

      const scal_region_info_t &info = region_[index];
      const item_t *sp_mem_reg = &data_array[info.sp_offset];
      const item_t *so_mem_reg = &data_array[info.so_offset];

      return scal_region_t (info,
                            data_vector_t (sp_mem_reg + placement_info_.sp_offset,   placement_info_.sp_step,    info.Sp_count),
                            data_vector_t (so_mem_reg + placement_info_.so_offset,   placement_info_.so_step,    info.Sp_count),
                            data_vector_t (sp_mem_reg + placement_info_.krp_offset,  placement_info_.krp_step,   info.So_count),
                            data_vector_t (so_mem_reg + placement_info_.krop_offset, placement_info_.krop_step,  info.So_count),
                            data_vector_t (sp_mem_reg + placement_info_.pcp_offset,  placement_info_.pcp_step,   info.Sp_count)
                           );
    }

    void
    init_regions ()
    {
      item_array_t &data_array = *data_;

      if (region_2_.size ())
        return ;

      for (size_t i = 0, cnt = region_.size (); i < cnt; ++i)
        {
          const scal_region_info_t &info  = region_[i];
          const item_t *sp_mem_reg        = &data_array[info.sp_offset];
          const item_t *so_mem_reg        = &data_array[info.so_offset];

        region_2_.push_back (new scal_region_t (info,
                              data_vector_t (sp_mem_reg + placement_info_.sp_offset,   placement_info_.sp_step,    info.Sp_count),
                              data_vector_t (so_mem_reg + placement_info_.so_offset,   placement_info_.so_step,    info.Sp_count),
                              data_vector_t (sp_mem_reg + placement_info_.krp_offset,  placement_info_.krp_step,   info.So_count),
                              data_vector_t (so_mem_reg + placement_info_.krop_offset, placement_info_.krop_step,  info.So_count),
                              data_vector_t (sp_mem_reg + placement_info_.pcp_offset,  placement_info_.pcp_step,   info.Sp_count)
                             ));
        }
    }

    const scal::data_placement::scal_placement_info &
    get_placement_info () const
    {
      return placement_info_;
    }

    const sp_array_item_t 
    get_raw_data () const
    {
      return data_;
    }

    const region_vector_t &
    get_region_info () const
    {
      return region_;
    }

    void
    update_so (const sp_scal_data_t &water_data);

  private:

    item_t 
    get_phase_sat (item_t spr) const
    {
      return 0.0;
    }
    item_t 
    get_oil_sat (item_t sorp)  const
    {
      return 0.0;
    }

    static void precalc (scal_region_info_t &info, const scal_region_t &region, scal::data_placement::scal_data_placement_type type, bool is_water);
    static void precalc_min_index_and_spr (const data_vector_t &main, const data_vector_t &slave, int &main_index, item_t &spr);
    static bool check (const scal_region_t &region, bool is_water);

  private:

    scal::data_placement::scal_placement_info    placement_info_;
    sp_array_item_t                       data_;
    region_vector_t                       region_;
    region_vector_2_t                     region_2_;

  public:

    BLUE_SKY_TYPE_DECL_T (scal_2p_data_holder);

  };

} 


#endif // #ifndef BS_SCAL_SCAL_2P_DATA_HOLDER_H_
