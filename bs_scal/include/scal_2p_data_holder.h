/**
 * \file scal_2p_data_holder.h
 * \brief holder for SCAL data
 * \author Sergey Miryanov
 * \date 30.12.2008
 * */
#ifndef BS_SCAL_SCAL_2P_DATA_HOLDER_H_
#define BS_SCAL_SCAL_2P_DATA_HOLDER_H_

#include "scal_data_placement_info.h"
#include "table_iface.h"

#include "bs_serialize_decl.h"

namespace blue_sky {

  /**
   * \brief hold scal region data (So, Sp, Krp, Krop)
   * \detail we hold data as a plain memory region and a list of region_info
   * */
  class BS_API_PLUGIN scal_2p_data_holder : public scal_2p_data_holder_iface
  {
  public:
    typedef t_float                                 item_t;
    typedef t_int                                   index_t;

    typedef scal_2p_data_holder                     this_t;

    typedef scal_region                             scal_region_t;
    typedef scal_region_info <item_t>               scal_region_info_t;
    typedef data_vector <item_t>                    data_vector_t;
    typedef v_float                                 item_array_t;
    typedef smart_ptr <item_array_t, true>          sp_array_item_t;
    typedef std::vector <scal_region_info_t>        region_vector_t;
    typedef std::vector <scal_region_t *>           region_vector_2_t;
    typedef smart_ptr <this_t, true>				sp_scal_data_t;
    typedef BS_SP (table_iface)                     sp_table;
    typedef std::vector <sp_table>                  sp_table_array_t;
    
    enum
      {
        SCAL_TABLE_SP = 0,
        SCAL_TABLE_SO,
        SCAL_TABLE_KRP,
        SCAL_TABLE_KROP,
        SCAL_TABLE_PCP,
        SCAL_TABLE_TOTAL
      };

	enum
	{
		SCAL_TABLE_SP2 = 0,
		SCAL_TABLE_KRP2,
		SCAL_TABLE_KROP2,
		SCAL_TABLE_PCP2
	};
    
    
    void add_spof (sp_array_item_t const &data, t_long region_index, bool is_water);
    void add_spfn (sp_array_item_t const &data, t_long region_index, bool is_water);
    void add_sof3 (sp_array_item_t const &data, t_long region_index, bool is_water);
    void add_sof2 (sp_array_item_t const &data, t_long region_index, bool is_water);

    ~scal_2p_data_holder ()
    {
      for (size_t i = 0, cnt = region_2_.size (); i < cnt; ++i)
        {
          delete region_2_[i];
        }
      scal_table_array.clear ();  
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
                            data_vector_t (so_mem_reg + placement_info_.so_offset,   placement_info_.so_step,    info.So_count),
                            data_vector_t (sp_mem_reg + placement_info_.krp_offset,  placement_info_.krp_step,   info.Sp_count),
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
                              data_vector_t (so_mem_reg + placement_info_.so_offset,   placement_info_.so_step,    info.So_count),
                              data_vector_t (sp_mem_reg + placement_info_.krp_offset,  placement_info_.krp_step,   info.Sp_count),
                              data_vector_t (so_mem_reg + placement_info_.krop_offset, placement_info_.krop_step,  info.So_count),
                              data_vector_t (sp_mem_reg + placement_info_.pcp_offset,  placement_info_.pcp_step,   info.Sp_count)
                             ));
        }
    }
 
	void
	clear_regions ()
	{
		region_.clear();
		region_2_.clear();
	}

    void 
    init_table_array (const t_long n_scal_regions)
    {
      BS_ASSERT (n_scal_regions >= 1);
      scal_table_array.resize (n_scal_regions);
      for (t_long i = 0; i < n_scal_regions; i++)
        {
          scal_table_array[i] = BS_KERNEL.create_object ("table");
		  scal_table_array[i]->init(0, SCAL_TABLE_TOTAL);
        }  
    }

	void init_regions_from_table(sp_table const &scal_table, bool is_water)
	{
		if (region_2_.size ())
			return ;
		scal_table_array.resize (1);
		scal_table_array[0] = scal_table;
		scal_region_info_t info;

        const t_int sp_table_len = scal_table->get_n_rows (SCAL_TABLE_SP);
        const t_int so_table_len = sp_table_len;  //scal_table->get_n_rows (SCAL_TABLE_SO);

		// do explicit conversion from doubles returned by scal_table
		// into t_float required by data_vector_t

		const t_double* sp__ = scal_table->get_col_ptr (SCAL_TABLE_SP2);
		const std::vector< t_float > vsp_(sp__, sp__ + sp_table_len);
		const item_t *sp_ = &vsp_[0];
		//const item_t *sp_   = scal_table->get_col_ptr (SCAL_TABLE_SP2);
          //const item_t *so_   = scal_table->get_col_ptr (SCAL_TABLE_SO); 

		const t_double* krp__ = scal_table->get_col_ptr (SCAL_TABLE_KRP2);
		const std::vector< t_float > vkrp_(krp__, krp__ + sp_table_len);
		const item_t *krp_  = &vkrp_[0];
		//const item_t *krp_  = scal_table->get_col_ptr (SCAL_TABLE_KRP2);

		const t_double* krop__ = scal_table->get_col_ptr (SCAL_TABLE_KROP2);
		const std::vector< t_float > vkrop_(krop__, krop__ + sp_table_len);
		const item_t *krop_  = &vkrop_[0];
		//const item_t *krop_ = scal_table->get_col_ptr (SCAL_TABLE_KROP2); 

		const t_double* pcp__ = scal_table->get_col_ptr (SCAL_TABLE_PCP2);
		const std::vector< t_float > vpcp_(pcp__, pcp__ + sp_table_len);
		const item_t *pcp_  = &vpcp_[0];
		//const item_t *pcp_  = scal_table->get_col_ptr (SCAL_TABLE_PCP2); 
          
        item_t *so_ = new item_t[so_table_len];
		item_t *_pcp_ = new item_t[sp_table_len];
		for (int i = 0; i < sp_table_len; i++)
		{
			so_[i] = 1.0 - sp_[i];
			_pcp_[i]  = is_water ? - pcp_[i] : pcp_[i];
		}

        region_2_.push_back (new scal_region_t (info,
                              data_vector_t (sp_,   1, sp_table_len),
                              data_vector_t (so_,   1, so_table_len),
                              data_vector_t (krp_,  1, sp_table_len),
                              data_vector_t (krop_, 1, so_table_len),
                              data_vector_t (_pcp_,  1, sp_table_len)
                             ));
	}

    void
    init_regions_from_tables ()
    {
      if (region_2_.size ())
        return ;
      
      BS_ASSERT (scal_table_array.size () > 0); 
      for (t_long i = 0, cnt = scal_table_array.size (); i < cnt; ++i)
        {
          const scal_region_info_t &info  = region_[i];
          sp_table scal_table = scal_table_array[i];
          
          const t_int sp_table_len = scal_table->get_n_rows (SCAL_TABLE_SP);
          const t_int so_table_len = scal_table->get_n_rows (SCAL_TABLE_SO);

          const t_double* sp__ = scal_table->get_col_ptr (SCAL_TABLE_SP);
          const std::vector< t_float > vsp_(sp__, sp__ + sp_table_len);
          const item_t *sp_ = &vsp_[0];

          const t_double* so__ = scal_table->get_col_ptr (SCAL_TABLE_SO);
          const std::vector< t_float > vso_(so__, so__ + so_table_len);
          const item_t *so_ = &vso_[0];

          const t_double* krp__ = scal_table->get_col_ptr (SCAL_TABLE_KRP);
          const std::vector< t_float > vkrp_(krp__, krp__ + sp_table_len);
          const item_t *krp_  = &vkrp_[0];
          
          const t_double* krop__ = scal_table->get_col_ptr (SCAL_TABLE_KROP);
          const std::vector< t_float > vkrop_(krop__, krop__ + so_table_len);
          const item_t *krop_  = &vkrop_[0];
          
          const t_double* pcp__ = scal_table->get_col_ptr (SCAL_TABLE_PCP);
          const std::vector< t_float > vpcp_(pcp__, pcp__ + sp_table_len);
          const item_t *pcp_  = &vpcp_[0];

          //const item_t *sp_   = scal_table->get_col_ptr (SCAL_TABLE_SP);
          //const item_t *so_   = scal_table->get_col_ptr (SCAL_TABLE_SO); 
          //const item_t *krp_  = scal_table->get_col_ptr (SCAL_TABLE_KRP); 
          //const item_t *krop_ = scal_table->get_col_ptr (SCAL_TABLE_KROP); 
          //const item_t *pcp_  = scal_table->get_col_ptr (SCAL_TABLE_PCP); 
          
          
          region_2_.push_back (new scal_region_t (info,
                              data_vector_t (sp_,   1, sp_table_len),
                              data_vector_t (so_,   1, so_table_len),
                              data_vector_t (krp_,  1, sp_table_len),
                              data_vector_t (krop_, 1, so_table_len),
                              data_vector_t (pcp_,  1, sp_table_len)
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
    get_phase_sat (item_t /*spr*/) const
    {
      return 0.0;
    }
    item_t 
    get_oil_sat (item_t /*sorp*/)  const
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
    
    sp_table_array_t                      scal_table_array;
  public:

    BLUE_SKY_TYPE_DECL_T (scal_2p_data_holder);

    friend class blue_sky::bs_serialize;

  };

} 


#endif // #ifndef BS_SCAL_SCAL_2P_DATA_HOLDER_H_
