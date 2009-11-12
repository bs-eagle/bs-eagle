/**
 *       \file  well_controller.h
 *      \brief  Well controller 
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  14.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
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

  /**
   * \class rate_data
   * \brief Stores rate data
   * */
  template <typename strategy_t>
  struct rate_data
  {
    typedef typename strategy_t::item_t item_t;

    /**
     * \class rate_data_inner
     * \brief Stores specific data for production 
     *        and injection rates
     * */
    struct rate_data_inner
    {
      item_t oil;
      item_t water;
      item_t gas;
      item_t liquid;

      //! ctor
      rate_data_inner ()
      : oil (0)
      , water (0)
      , gas (0)
      , liquid (0)
      {
      }

      //! Sets data to value
      void
      operator= (item_t value)
      {
        oil     = value;
        water   = value;
        gas     = value;
        liquid  = value;
      }

      //! Sum two rate_data_inner objects
      void
      operator+= (const rate_data_inner &rhs)
      {
        oil     += rhs.oil;
        water   += rhs.water;
        gas     += rhs.gas;
        liquid  += rhs.liquid;
      }

      /**
       * \brief  Multiplies data on mult
       * \param  mult
       * \return New rate_data_inner
       * */
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

    //! ctor
    rate_data ()
    : liquid_inner (0)
    , free_gas (0)
    , solution_gas (0)
    {
    }

    /**
     * \brief  Sets data to value
     * \param  value
     * \return Reference to this object
     * */
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

    /**
     * \brief  Sum two rate_data objects
     * \param  rhs
     * */
    void
    operator += (const rate_data <strategy_t> &rhs)
    {
      prod          += rhs.prod;
      inj           += rhs.inj;
      liquid_inner  += rhs.liquid_inner;
      free_gas      += rhs.free_gas;
      solution_gas  += rhs.solution_gas;
    }

    /**
     * \brief  Multiplies data on mult
     * \param  mult
     * \return New rate_data
     * */
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

    rate_data_inner prod; //!< Production part of rate data
    rate_data_inner inj;  //!< Injection part of rate data

    item_t liquid_inner;
    item_t free_gas;      //!< Free gas
    item_t solution_gas;  //!< Solution (?) gas
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
    /**
     * \enum  rate_value_type
     * \brief Type of rate value
     * */
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

    /**
     * \brief  Returns true if type is oil_rate_value
     * \param  type
     * \return True if type is oil_rate_value
     * */
    bool 
    is_oil_rate_value (rate_value_type type);

    /**
     * \brief  Returns true if type is water_rate_value
     * \param  type
     * \return True if type is water_rate_value
     * */
    bool 
    is_water_rate_value (rate_value_type type);

    /**
     * \brief  Returns true if type is gas_rate_value
     * \param  type
     * \return True if type is gas_rate_value
     * */
    bool 
    is_gas_rate_value (rate_value_type type);

    /**
     * \brief  Returns true if type is oil_rate_value
     *         or is water_rate_value
     * \param  type
     * \return True if type is oil_rate_value
     *         or is water_rate_value
     * */
    bool 
    is_liquid_rate_value (rate_value_type type);

    /**
     * \brief  Converts string value to rate_value_type
     * \param  str String value to convert
     * \return Throws exception if value is invalid
     *         otherwise element of rate_value_type
     * */
    rate_value_type
    rate_value_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    /**
     * \enum  injection_type
     * \brief Type of injection
     * */
    enum injection_type
    {
      injection_none,     //!< Invalid injection type
      injection_water,    //!< Water injection
      injection_gas,      //!< Gas injection
      injection_oil,      //!< Oil injection
    };

    /**
     * \brief  Converts string value to injection_type
     * \param  str String value to convert
     * \return Throws exception if value is invalid
     *         otherwise element of injection_type
     * */
    injection_type
    injection_type_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    /**
     * \class well_controller
     * \brief Well controller
     * */
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
        //! dtor
        virtual ~well_controller () {}

        //! Clears rate
        void 
        clear_rate ();

        /**
         * \brief  Sets rate value
         * \param  rate_value Type of rate to set
         * \param  value Value to set
         * */
        void 
        set_rate (rate_value_type rate_value, item_t value);

        /**
         * \brief  Sets BHP value
         * \param  value BHP value
         * */
        void 
        set_bhp (item_t value);

        /**
         * \brief  Sets historical value of BHP
         * \param  value Value to set
         * */
        void 
        set_bhp_history (item_t value);

        /**
         * \brief  Adds BHP control (controls well by BHP)
         * \param  bhp_control 
         * */
        void 
        add_bhp_control (const sp_rate_control_t &bhp_control);

        /**
         * \brief  Adds rate control (controls well by rate)
         * \param  rate_control
         * */
        void 
        add_rate_control (const sp_rate_control_t &rate_control);

        /**
         * \brief  Sets main control (bhp or rate depends on control value)
         * \param  well
         * \param  control
         * */
        void 
        set_main_control (const sp_well_t &well, rate_control_type control);

        /**
         * \brief  Sets injection type
         * \param  type
         * */
        void 
        set_injection_type (injection_type type);

        /**
         * \brief  Returns rate data
         * \return Rate data
         * */
        const rate_data_t &
        rate () const;

        /**
         * \brief  Returns BHP value
         * \return BHP value
         * */
        const item_t &
        bhp () const;

        /**
         * \brief  Returns historical BHP value
         * \return Historical BHP
         * */
        const item_t &
        bhp_history () const;

        /**
         * \brief  Returns injection type
         * \return Injection type
         * */
        const injection_type &
        injection () const;

        /**
         * \brief  Checks is control is bhp
         * \return True if control is bhp_control
         * */
        bool 
        is_bhp () const;

        /**
         * \brief  Checks is control is rate
         * \return True if control is rate_control
         * */
        bool 
        is_rate () const;

        /**
         * \brief  Checks is well working in production mode
         * \return True if well is production
         * */
        bool 
        is_production () const;

        /**
         * \brief  Checks is given BHP value is valid
         * \param  pressure Value to be compared with BHP value
         * \param  bhp Value to be checked
         * \return True if value if valid
         * */
        bool 
        is_valid_connection_bhp (item_t pressure, item_t bhp) const;

        /**
         * \brief  Saves well_controller internal state
         * */
        void 
        save_control ();
        /**
         * \brief  Restores well_controller internal state
         * \return True is restored successfully
         * */
        bool 
        restore_control ();
        /**
         * \brief  Saves well_controller internal state
         * */
        void 
        save_niter_control ();
        /**
         * \brief  Restores well_controller internal state
         * \return True is restored successfully
         * */
        bool 
        restore_niter_control ();

        /**
         * \brief  Switches well to control by BHP
         * \param  well Well to be swithced
         * */
        void 
        switch_to_bhp (sp_well_t &well);
        /**
         * \brief  Checks is well works properly
         * \param  well Well to be checked
         * \return True if control of well was changed due checks
         * */
        bool 
        check (sp_well_t &well);

        /**
         * \brief  Calculates rate
         * \todo   Obsolete, should be removed
         * */
        void 
        calc_rate (const sp_calc_model_t &calc_model, sp_well_t &well, sp_jmatrix_t &jmatrix) const;
        /**
         * \brief  Calculates derivs
         * \todo   Obsolete, should be removed
         * */
        void 
        calc_derivs (const sp_calc_model_t &calc_model, sp_well_t &well, sp_jmatrix_t &jmatrix) const;

        /**
         * \brief  Returns type of control
         * \return Type of control
         * */
        rate_control_type 
        get_control_type () const;

        //! blue-sky type declaration
        BLUE_SKY_TYPE_DECL_T (well_controller);

      public:
        rate_data_t                     rate_;                //!< Rates

      private:
        item_t                          bhp_;                 //!< pw,ref; rate_value_type::pref_value
        item_t                          bhp_history_;         //!<

        injection_type                  injection_type_;      //!< Injection type (now only WATER injection supports)

        static sp_rate_control_t        dummy_control_;
        sp_rate_control_t               bhp_control_;
        sp_rate_control_t               rate_control_;

        sp_rate_control_t               current_control_;
        sp_rate_control_t               saved_control_;
        sp_rate_control_t               saved_niter_control_;
      };

    /**
     * \class well_controller_factory
     * \brief Factory of well_controllers
     * \todo  Obsolete, should be redisigned
     * */
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

        //! dtor
        virtual ~well_controller_factory () {};

        /**
         * \brief  Sets pointer to factory of rate_control
         * \param  rate_control_factory
         * */
        void 
        set_rate_control_factory (const sp_well_rate_control_factory_t &rate_control_factory);

        /**
         * \brief  Creates well_controller
         * \return Instance of well_controller
         * */
        virtual sp_well_controller_t
        create_controller () const;

        /**
         * \brief  Creates well_controller_factory
         * \param  rate_control
         * \param  is_prod
         * \param  calc_model
         * \return Instance of rate_control
         * */
        virtual sp_rate_control_t
        create_control (rate_control_type rate_control, bool is_prod, const sp_calc_model_t &calc_model) const;

        //! blue-sky type declaration
        BLUE_SKY_TYPE_DECL_T (well_controller_factory);

      private:

        sp_well_rate_control_factory_t well_rate_control_factory_; //!< rate_control factory
      };

    /**
     * \brief  Registers well_controller types in blue-sky kernel
     * \param  pd plugin_descriptor
     * \return True if all types registered successfully
     * */
    bool
    well_controller_register_type (const blue_sky::plugin_descriptor &pd);

    /**
     * \brief  Registers well_controller_factory types in blue-sky kernel
     * \param  pd plugin_descriptor
     * \return True if all types registered successfully
     * */
    bool
    well_controller_factory_register_type (const blue_sky::plugin_descriptor &pd);


  } // namespace wells
} // namespace blue_sky


#endif  // #ifnded BS_WELL_CONTROLLER_H_

