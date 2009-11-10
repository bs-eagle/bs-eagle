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
  template <typename strategy_t>
  class BS_API_PLUGIN scal_2p_data_holder : public objbase
  {
  public:
    typedef typename strategy_t::item_t             item_t;
    typedef typename strategy_t::index_t            index_t;

    typedef scal_2p_data_holder<strategy_t>         this_t;

    typedef scal_region <strategy_t>                scal_region_t;
    typedef scal_region_info <strategy_t>           scal_region_info_t;
    typedef data_vector<strategy_t>                 data_vector_t;
    typedef typename strategy_t::item_array_t       item_array_t;
    typedef seq_vector <scal_region_info_t>         region_vector_t;
    typedef seq_vector <scal_region_t *>            region_vector_2_t;

    typedef smart_ptr <this_t, true>								sp_scal_data_t;

    void add_spof (const item_array_t &data, bool is_water);
    void add_spfn (const item_array_t &data, size_t region_index, bool is_water);
    void add_sof3 (const item_array_t &data, size_t region_index, bool is_water);

    ~scal_2p_data_holder ()
    {
      for (size_t i = 0, cnt = region_2_.size (); i < cnt; ++i)
        {
          delete region_2_[i];
        }
    }

    const scal_region_t &
    get_region (index_t index) const
    {
      return *region_2_[index];
    }
    scal_region_t 
    get_region_internal (int index) const
    {
      BS_ASSERT (index >= 0 && index < (int)region_.size ());

      const scal_region_info_t &info = region_[index];
      const item_t *sp_mem_reg = &data_[info.sp_offset];
      const item_t *so_mem_reg = &data_[info.so_offset];

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
      if (region_2_.size ())
        return ;

      for (size_t i = 0, cnt = region_.size (); i < cnt; ++i)
        {
          const scal_region_info_t &info  = region_[i];
          const item_t *sp_mem_reg        = &data_[info.sp_offset];
          const item_t *so_mem_reg        = &data_[info.so_offset];

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

    const item_array_t &
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

    scal::data_placement::scal_placement_info
    placement_info_;
    item_array_t                          data_;
    region_vector_t                       region_;
    region_vector_2_t                     region_2_;

  public:

    BLUE_SKY_TYPE_DECL_T (scal_2p_data_holder);

  };

} 


#endif // #ifndef BS_SCAL_SCAL_2P_DATA_HOLDER_H_
