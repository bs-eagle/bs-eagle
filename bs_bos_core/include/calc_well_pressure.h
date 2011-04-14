/**
 *       \file  calc_well_pressure.h
 *      \brief  Base class for objects that calculate BHP for well
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  26.09.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_CALC_WELL_PRESSURE_H_
#define BS_CALC_WELL_PRESSURE_H_

//// WTF??
//#include "well_results_storage.h"
//#include "fip_results_storage.h"

namespace blue_sky
  {

  class calc_model;
  class well;

  /**
   * \class calc_well_pressure_base
   * \brief Base class for objects that calculate BHP for well
   * \todo  Should be renamed to calc_well_pressure_iface
   * */
  class calc_well_pressure_base : public objbase
    {
    public:
      typedef t_long                            index_t;
      typedef t_double                          item_t;

      typedef calc_model                        calc_model_t;
      typedef well                              well_t;
      typedef smart_ptr <calc_model_t, true>    sp_calc_model_t;
      typedef smart_ptr <well_t, true>          sp_well_t;

    public:
      /**
       * \brief  calc_well_pressure_base dtor
       * \param  
       * \return 
       * */
      virtual ~calc_well_pressure_base () {}

      /**
       * \brief  Calculates BHP for well
       * \param  well
       * \param  calc_model
       * \return True if well should be switched to controll by BHP otherwise false
       * */
      virtual bool 
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const = 0;
    };

  /**
   * \class calc_well_pressure
   * \brief Calculates BHP for well
   * */
  class calc_well_pressure : public calc_well_pressure_base 
    {
    public:
      typedef t_long                    index_t;
      typedef t_double                  item_t;
      typedef v_double                  item_array_t;

      typedef calc_well_pressure_base   base_t;
      typedef calc_well_pressure        this_t;

      typedef base_t::sp_calc_model_t   sp_calc_model_t;
      typedef base_t::sp_well_t         sp_well_t;

    public:

      /**
       * \brief  Calculates BHP for well
       * \param  well
       * \param  calc_model
       * \return True if well should be switched to controll by BHP otherwise false
       * */
      bool
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const;

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL (calc_well_pressure);

    protected:

      /**
       * \brief  Calculates BHP for well if well controlled by rate
       * \param  well
       * \param  calc_model
       * \return True if well should be switched to controll by BHP otherwise false
       * */
      bool
      calculate_for_rate (sp_well_t &well, const sp_calc_model_t &calc_model) const;

    };

  /**
   * \brief  Registers calc_well_pressure types in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return True if all types registered successfully
   * */
  bool
  calc_well_pressure_register_types (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky

#endif  // #ifndef BS_CALC_WELL_PRESSURE_H_

