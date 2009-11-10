#ifndef PY_CALC_MODEL_H
#define PY_CALC_MODEL_H

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_pvt.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  struct calc_model_data;

  namespace python
    {

    template <typename strategy_t>
    struct py_calc_model_data : calc_model_data <strategy_t>
      {
        typedef py_calc_model_data <strategy_t>       this_t;
        typedef calc_model_data <strategy_t>          base_t;
        typedef typename strategy_t::item_t           item_t;
        typedef typename strategy_t::item_array_t     item_array_t;
        typedef typename strategy_t::index_t          index_t;

        py_calc_model_data () : base_t () {}
        py_calc_model_data (const base_t &cm) : base_t (cm) {}

        item_array_t get_cap_pressure () const;
        item_array_t get_s_deriv_cap_pressure () const;

        item_array_t get_relative_perm () const;
        item_array_t get_s_deriv_relative_perm () const;

        //item_t get_gas_oil_ratio () const;
        item_t get_p_deriv_gas_oil_ratio () const;

        item_array_t get_invers_fvf () const;
        item_array_t get_p_deriv_invers_fvf () const;
        item_t get_gor_deriv_invers_fvf () const;

        item_array_t get_invers_viscosity () const;
        item_array_t get_p_deriv_invers_viscosity () const;
        item_t get_gor_deriv_invers_viscosity () const;

        item_array_t get_invers_visc_fvf () const;
        item_array_t get_p_deriv_invers_visc_fvf () const;
        item_t get_gor_deriv_invers_visc_fvf () const;

        item_array_t get_density () const;
        item_array_t get_p_deriv_density () const;
        item_t get_gor_deriv_density () const;

        item_t get_porosity () const;
        item_t get_p_deriv_porosity () const;

        item_t get_truns_mult () const;
        item_t get_p_deriv_truns_mult () const;

        item_array_t get_mobility () const;
        item_array_t get_p_deriv_mobility () const;
        item_array_t get_s_deriv_mobility () const;

        item_array_t get_prev_fluid_volume () const;
      };

    template <typename strategy_t>
    class py_scal_3p;

    template <typename strategy_t>
    class py_idata;

    template<class strategy_t>
    class py_calc_model : public py_objbase
      {
      public:
        typedef typename strategy_t::item_array_t   item_array_t;
        typedef typename strategy_t::item_t         item_t;
        typedef typename strategy_t::index_array_t  index_array_t;
        typedef typename strategy_t::index_t        index_t;

        typedef calc_model<strategy_t>              wrapped_t;
        typedef smart_ptr<wrapped_t>					      sp_cm_t;

        typedef py_idata <strategy_t>               py_idata_t;

        typedef py_calc_model_data <strategy_t>     py_calc_model_data_t;

        typedef typename wrapped_t::sp_fi_params    sp_fi_params_t;

      public:

        py_calc_model ();
        py_calc_model(const sp_cm_t &);
        ~py_calc_model () {}

        py_pvt_water <pvt_water <strategy_t> >
        get_pvt_water (index_t n_pvt_region) const;

        py_pvt_gas <pvt_gas <strategy_t> >
        get_pvt_gas (index_t n_pvt_region) const;

        py_pvt_oil <pvt_oil <strategy_t> >
        get_pvt_oil (index_t n_pvt_region) const;

        // methods
        item_array_t get_saturation ();
        item_array_t get_pressure();
        item_array_t get_gas_oil_ratio();
        index_t get_pvt_num();
        index_t get_sat_num();
        index_t get_fip_num();

        item_array_t get_plane_flow_rate ();
        item_array_t get_full_step_plane_flow_rate ();

        index_array_t get_pvt_regions ();
        index_array_t get_sat_regions ();
        index_array_t get_fip_regions ();
        index_array_t get_rock_regions ();

        index_array_t get_bconn_mainvar ();

        item_array_t get_bconn_pressure ();
        item_array_t get_bconn_saturation ();
        item_array_t get_bconn_gor ();

        py_scal_3p<strategy_t>
        get_scal () const;

        py_calc_model_data_t get_data (index_t) const;

        index_t cmd_num () const;
        index_t all_data_len () const;

        index_t cap_pressure_len () const;
        item_t get_cap_pressure (index_t, index_t) const;

        index_t s_deriv_cap_pressure_len () const;
        item_t get_s_deriv_cap_pressure (index_t, index_t) const;

        index_t relative_perm_len () const;
        item_t get_relative_perm (index_t, index_t) const;

        index_t s_deriv_relative_perm_len () const;
        item_t get_s_deriv_relative_perm (index_t, index_t) const;

        index_t p_deriv_gas_oil_ratio_len () const;
        item_t get_p_deriv_gas_oil_ratio (index_t, index_t) const;

        index_t invers_fvf_len () const;
        item_t get_invers_fvf (index_t, index_t) const;

        index_t p_deriv_invers_fvf_len () const;
        item_t get_p_deriv_invers_fvf (index_t, index_t) const;

        index_t gor_deriv_invers_fvf_len () const;
        item_t get_gor_deriv_invers_fvf (index_t, index_t) const;

        index_t invers_viscosity_len () const;
        item_t get_invers_viscosity (index_t, index_t) const;

        index_t p_deriv_invers_viscosity_len () const;
        item_t get_p_deriv_invers_viscosity (index_t, index_t) const;

        index_t gor_deriv_invers_viscosity_len () const;
        item_t get_gor_deriv_invers_viscosity (index_t, index_t) const;

        index_t invers_visc_fvf_len () const;
        item_t get_invers_visc_fvf (index_t, index_t) const;

        index_t p_deriv_invers_visc_fvf_len () const;
        item_t get_p_deriv_invers_visc_fvf (index_t, index_t) const;

        index_t gor_deriv_invers_visc_fvf_len () const;
        item_t get_gor_deriv_invers_visc_fvf (index_t, index_t) const;

        index_t density_len () const;
        item_t get_density (index_t, index_t) const;

        index_t p_deriv_density_len () const;
        item_t get_p_deriv_density (index_t, index_t) const;

        index_t gor_deriv_density_len () const;
        item_t get_gor_deriv_density (index_t, index_t) const;

        index_t porosity_len () const;
        item_t get_porosity (index_t, index_t) const;

        index_t p_deriv_porosity_len () const;
        item_t get_p_deriv_porosity (index_t, index_t) const;

        index_t truns_mult_len () const;
        item_t get_truns_mult (index_t, index_t) const;

        index_t p_deriv_truns_mult_len () const;
        item_t get_p_deriv_truns_mult (index_t, index_t) const;

        index_t mobility_len () const;
        item_t get_mobility (index_t, index_t) const;

        index_t p_deriv_mobility_len () const;
        item_t get_p_deriv_mobility (index_t, index_t) const;

        index_t s_deriv_mobility_len () const;
        item_t get_s_deriv_mobility (index_t, index_t) const;

        index_t prev_fluid_volume_len () const;
        item_t get_prev_fluid_volume (index_t, index_t) const;

        void initialize_datas ();

        int get_n_phases () const;
        bool is_water () const;
        bool is_gas () const;
        bool is_oil () const;

        sp_fi_params_t get_fi_params () const;

      private:
      };

    void py_export_calc_model ();
  }
}

#endif // PY_CALC_MODEL_H
