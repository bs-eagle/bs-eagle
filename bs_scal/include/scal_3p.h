/*
 * \file scal_3p.h
 * \brief scal 3p
 * \author Sergey Miryanov
 * \date 19.05.2008
 * */
#ifndef BS_SCAL_3P_H_
#define BS_SCAL_3P_H_

#include "calc_model_data.h"

namespace blue_sky
{

  template <typename strategy_t>
  class BS_API_PLUGIN scal_region;

  template <typename strategy_t>
  class BS_API_PLUGIN scale_array_holder;

  template <typename strategy_t>
  class BS_API_PLUGIN scal_2p_data_holder;

  template <typename strategy_t>
  class BS_API_PLUGIN jfunction;

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  class BS_API_PLUGIN scal_3p : public objbase
    {
    public:
      typedef scal_3p <strategy_t>											this_t;

      typedef typename strategy_t::item_t								item_t;
      typedef typename strategy_t::index_t							index_t;
      typedef typename strategy_t::index_array_t        index_array_t;
      typedef typename strategy_t::item_array_t         item_array_t;

      typedef scal_region <strategy_t>									scal_region_t;
      typedef scale_array_holder <strategy_t>						scale_array_holder_t;
      typedef scal_2p_data_holder <strategy_t>					scal_2p_data_holder_t;

      typedef jfunction <strategy_t>										jfunction_t;

      typedef smart_ptr <jfunction_t, true>							sp_jfunction_t;
      typedef smart_ptr <scale_array_holder_t, true>		sp_scale_array_holder_t;
      typedef smart_ptr <scal_2p_data_holder_t, true>		sp_scal_2p_data_holder_t;
      typedef boost::array <index_t, FI_PHASE_TOT>			phase_d_t;
      typedef boost::array <index_t, FI_PHASE_TOT>			sat_d_t;

      typedef calc_model_data <strategy_t>              data_t;
      typedef shared_vector <data_t>                    data_array_t;

      typedef unsigned char															phase_index_t;
      typedef unsigned char															sat_index_t;

      struct scal_3p_impl_base;

    public:

      ~scal_3p ();

      sp_scale_array_holder_t   
      get_water_scale () const
      {
        return water_scale;
      }
      sp_scale_array_holder_t
      get_gas_scale () const
      {
        return gas_scale;
      }

      sp_scal_2p_data_holder_t
      get_water_data () const
      {
        return water_data;
      }
      sp_scal_2p_data_holder_t
      get_gas_data () const
      {
        return gas_data;
      }

      sp_jfunction_t
      get_water_jfunction () const
      {
        return water_jfunc;
      }
      sp_jfunction_t
      get_gas_jfunction () const
      {
        return gas_jfunc;
      }

      void
      get_relative_perm (index_t cell_index, 
        const item_array_t &saturation, 
        const index_array_t &sat_regions, 
        item_array_t &relative_perm, 
        item_array_t &s_deriv_relative_perm) const;

        void
        get_capillary (index_t cell_index, 
          const item_array_t &saturation, 
          const index_array_t &sat_regions, 
          const item_array_t &perm, 
          const item_array_t &poro, 
          item_array_t &cap, 
          item_array_t &s_deriv_cap) const;

      void
      process (const item_array_t &saturation, 
        const index_array_t &sat_regions,
        const item_array_t &perm,
        const item_array_t &poro,
        data_array_t &data) const;

      void
      process_init (index_t cell_index, 
        const item_t *pressure, 
        index_t sat_reg, 
        const item_t *perm_array, 
        item_t poro,
        item_t *sat, 
        item_t *pc_limit) const;

      void
      calc_pcp (index_t cell_index, 
        const item_t sat, 
        index_t sat_reg, 
        item_t cap, 
        item_t &pcp) const;

      void
      calc_gas_water_zone (index_t cell_index, 
        index_t sat_reg, 
        const item_t *perm_array, 
        item_t poro, 
        item_t pcgw,
        item_t &sw, 
        item_t &sg) const;

      void
      init (bool is_w, bool is_g, bool is_o, 
        const phase_d_t &phase_d, const phase_d_t &sat_d,
        RPO_MODEL_ENUM rpo_model);

      void
      set_water_jfunction (sp_jfunction_t jfunc);

      void
      set_gas_jfunction (sp_jfunction_t jfunc);

      void
      update_gas_data ();

    private:

      sp_scal_2p_data_holder_t  water_data;
      sp_scal_2p_data_holder_t  gas_data;

      sp_scale_array_holder_t   water_scale;
      sp_scale_array_holder_t   gas_scale;

      sp_jfunction_t            water_jfunc;
      sp_jfunction_t            gas_jfunc;

      scal_3p_impl_base         *impl_;

    public:

      BLUE_SKY_TYPE_DECL_T (scal_3p);
    };

  bool scal_register_types (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky

#endif	// #ifndef BS_SCAL_3P_H_
