/**
 * \file well_controller.h
 * \brief well_controller
 * \author Sergey Miryanov
 * \date 14.07.2008
 * */
#ifndef BS_WELL_CONTROLLER_H_
#define BS_WELL_CONTROLLER_H_

#include "well_type_helper.h"
#include "rate_control_type.h"
#include "well_rate_control_interface.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class well;

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  struct calc_model_data;

  template <typename strategy_t>
  struct rate_data
  {
    typedef typename strategy_t::item_t item_t;

    struct rate_data_inner
    {
      item_t oil;
      item_t water;
      item_t gas;
      item_t liquid;

      rate_data_inner ()
      : oil (0)
      , water (0)
      , gas (0)
      , liquid (0)
      {
      }

      void
      operator= (item_t value)
      {
        oil     = value;
        water   = value;
        gas     = value;
        liquid  = value;
      }

      void
      operator+= (const rate_data_inner &rhs)
      {
        oil     += rhs.oil;
        water   += rhs.water;
        gas     += rhs.gas;
        liquid  += rhs.liquid;
      }

      rate_data_inner
      operator * (item_t mult) const
      {
        rate_data_inner r;
        r.oil     = mult * oil;
        r.water   = mult * water;
        r.gas     = mult * gas;
        r.liquid  = mult * liquid;

        return r;
      }
    };

    rate_data ()
    : liquid_inner (0)
    , free_gas (0)
    , solution_gas (0)
    {
    }

    rate_data <strategy_t> &
    operator= (item_t value)
    {
      prod          = value;
      inj           = value;
      liquid_inner  = value;
      free_gas      = value;
      solution_gas  = value;

      return *this;
    }

    void
    operator += (const rate_data <strategy_t> &rhs)
    {
      prod          += rhs.prod;
      inj           += rhs.inj;
      liquid_inner  += rhs.liquid_inner;
      free_gas      += rhs.free_gas;
      solution_gas  += rhs.solution_gas;
    }

    rate_data 
    operator * (item_t mult) const
    {
      rate_data r;
      r.prod          = prod * mult;
      r.inj           = inj * mult;
      r.liquid_inner  = liquid_inner * mult;
      r.free_gas      = free_gas * mult;
      r.solution_gas  = solution_gas * mult;

      return r;
    }

    rate_data_inner prod;
    rate_data_inner inj;

    item_t liquid_inner;
    item_t free_gas;
    item_t solution_gas;
  };

  namespace wells
    {

    template <typename strategy_t>
    class well_controller;

    template <typename strategy_t>
    class connection;

    template <typename strategy_t>
    class well_rate_control;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    enum rate_value_type
    {
      null_value              = 0,          //!< null value; means value not set
      oil_rate_value          = 1,          //!< value for oil rate
      water_rate_value        = 2,          //!< value for water rate
      liquid_rate_value       = 3,          //!< value for liquid rate
      gas_rate_value          = 4,          //!< value for gas rate
      rate_value              = 5,          //!< rate value for injection well (actual rate selected by injection phase)
      bhp_value               = 8,          //!< value of pressure on surface
      liquid_inner_rate_value = 9,          //!< value for liquid rate into bed
    };
    bool is_oil_rate_value (rate_value_type type);
    bool is_water_rate_value (rate_value_type type);
    bool is_gas_rate_value (rate_value_type type);
    bool is_liquid_rate_value (rate_value_type type);

    rate_value_type
    rate_value_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    enum injection_type
    {
      injection_none,
      injection_water,
      injection_gas,
      injection_oil,
    };

    injection_type
    injection_type_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    template <typename strategy_t>
    class BS_API_PLUGIN well_controller : public objbase
      {
      public:

        typedef typename strategy_t::item_t               item_t;
        typedef typename strategy_t::index_t              index_t;
        typedef typename strategy_t::item_array_t         item_array_t;
        typedef typename strategy_t::index_array_t        index_array_t;
        typedef well_rate_control <strategy_t>            well_rate_control_t;
        typedef well <strategy_t>                         well_t;
        typedef connection <strategy_t>                   connection_t;
        typedef jacobian_matrix <strategy_t>              jacobian_matrix_t;
        typedef well_controller <strategy_t>              this_t;
        typedef calc_model <strategy_t>                   calc_model_t;

        typedef rate_data <strategy_t>                    rate_data_t;

        typedef smart_ptr <calc_model_t, true>            sp_calc_model_t;
        typedef smart_ptr <jacobian_matrix_t, true>       sp_jmatrix_t;
        typedef smart_ptr <well_t, true>                  sp_well_t;

        typedef smart_ptr <well_rate_control_t, true>     sp_rate_control_t;
        typedef smart_ptr <connection_t, true>            sp_connection_t;
        typedef smart_ptr <jacobian_matrix_t, true>       sp_jacobian_matrix_t;
        typedef smart_ptr <this_t, true>                  sp_this_t;

        typedef wells::type_helper <strategy_t>           helper_t;

        typedef typename helper_t::item_rhs_block_t       item_rhs_block_t;
        typedef typename helper_t::item_ww_block_t        item_ww_block_t;
        typedef typename helper_t::item_q_rate_t          item_q_rate_t;
        typedef typename helper_t::item_q_rate_inflow_t   item_q_rate_inflow_t;
        typedef typename helper_t::item_gas_rate_t        item_gas_rate_t;

      public:
        virtual ~well_controller () {}

        void clear_rate ();
        void set_rate (rate_value_type rate_value, item_t value);
        void set_bhp (item_t value);
        void set_bhp_history (item_t value);
        void add_bhp_control (const sp_rate_control_t &bhp_control);
        void add_rate_control (const sp_rate_control_t &rate_control);
        void set_main_control (const sp_well_t &well, rate_control_type control);
        void set_injection_type (injection_type type);

        const rate_data_t &rate () const;
        const item_t &bhp () const;
        const item_t &bhp_history () const;

        const injection_type &injection () const;

        bool is_bhp () const;
        bool is_rate () const;
        bool is_production () const;
        bool is_valid_connection_bhp (item_t pressure, item_t bhp) const;

        void save_control ();
        bool restore_control ();
        void save_niter_control ();
        bool restore_niter_control ();

        void switch_to_bhp (sp_well_t &well);
        bool check (sp_well_t &well);
        void calc_rate (const sp_calc_model_t &calc_model, sp_well_t &well, sp_jmatrix_t &jmatrix) const;
        void calc_derivs (const sp_calc_model_t &calc_model, sp_well_t &well, sp_jmatrix_t &jmatrix) const;

        rate_control_type get_control_type () const;

        BLUE_SKY_TYPE_DECL_T (well_controller);

      public:
        // TODO:
        rate_data_t                     rate_;

      private:
        item_t                          bhp_;                 //!< pw,ref; rate_value_type::pref_value
        item_t                          bhp_history_;         //!<

        injection_type                  injection_type_;      //!< injection type (now only WATER)

        static sp_rate_control_t        dummy_control_;
        sp_rate_control_t               bhp_control_;
        sp_rate_control_t               rate_control_;

        sp_rate_control_t               current_control_;
        sp_rate_control_t               saved_control_;
        sp_rate_control_t               saved_niter_control_;
      };

    template <typename strategy_t>
    class BS_API_PLUGIN well_controller_factory : public objbase
      {
      public:

        typedef calc_model <strategy_t>                       calc_model_t;
        typedef well_controller <strategy_t>                  well_controller_t;
        typedef well_rate_control <strategy_t>                well_rate_control_t;
        typedef well_rate_control_factory <strategy_t>        well_rate_control_factory_t;

        typedef smart_ptr <calc_model_t, true>                sp_calc_model_t;
        typedef smart_ptr <well_controller_t, true>           sp_well_controller_t;
        typedef smart_ptr <well_rate_control_t, true>         sp_rate_control_t;
        typedef smart_ptr <well_rate_control_factory_t, true> sp_well_rate_control_factory_t;

      public:

        virtual ~well_controller_factory () {};

        void set_rate_control_factory (const sp_well_rate_control_factory_t &rate_control_factory);

        virtual sp_well_controller_t    create_controller () const;
        virtual sp_rate_control_t           create_control (rate_control_type rate_control, bool is_prod, const sp_calc_model_t &calc_model) const;

        BLUE_SKY_TYPE_DECL_T (well_controller_factory);

      private:

        sp_well_rate_control_factory_t well_rate_control_factory_;
      };

    bool
    well_controller_register_type (const blue_sky::plugin_descriptor &pd);

    bool
    well_controller_factory_register_type (const blue_sky::plugin_descriptor &pd);


  } // namespace wells
} // namespace blue_sky


#endif  // #ifnded BS_WELL_CONTROLLER_H_

