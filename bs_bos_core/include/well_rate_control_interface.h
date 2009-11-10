/**
 * \file well_rate_control_interface.h
 * \brief interface for well_rate_control_impl
 * \author Sergey Miryanov
 * \date 21.11.2008
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_INTERFACE_H_
#define BS_WELLS_WELL_RATE_CONTROL_INTERFACE_H_

#include "rate_control_type.h"
#include "make_me_happy.h"

namespace blue_sky {

  template <typename strategy_t> class calc_model;
  template <typename strategy_t> class BS_API_PLUGIN well;

  namespace wells {

    template <typename strategy_t> class BS_API_PLUGIN connection;
    template <typename strategy_t> class BS_API_PLUGIN well_controller;
  }

}

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  class well_rate_control_interface : public objbase
  {
  public:
    typedef typename strategy_t::item_t         item_t;
    typedef typename strategy_t::index_t        index_t;

    typedef calc_model <strategy_t>             calc_model_t;
    typedef jacobian_matrix <strategy_t>        jmatrix_t;
    typedef well <strategy_t>                   well_t;
    typedef wells::well_controller <strategy_t> well_controller_t;
    typedef wells::connection <strategy_t>      connection_t;

    typedef smart_ptr <calc_model_t, true>      sp_calc_model_t;
    typedef smart_ptr <jmatrix_t, true>         sp_jmatrix_t;
    typedef smart_ptr <well_t, true>            sp_well_t;
    typedef smart_ptr <well_controller_t, true> sp_well_controller_t;
    typedef smart_ptr <connection_t, true>      sp_connection_t;

  public:
    virtual void compute_rate (const sp_calc_model_t &calc_model, sp_jmatrix_t &jmatrix, sp_well_t &well, const sp_well_controller_t &well_controller) const = 0;
    virtual void compute_derivs (const sp_calc_model_t &calc_model, sp_jmatrix_t &jmatrix, sp_well_t &well, const sp_well_controller_t &well_controller) const = 0;

    virtual ~well_rate_control_interface ()
    {
    }
  };

  template <typename strategy_t>
  class well_rate_control_factory : public objbase
  {
  public:
    typedef calc_model <strategy_t>                   calc_model_t;
    typedef well_rate_control_interface <strategy_t>  well_rate_control_t;

    typedef smart_ptr <calc_model_t, true>            sp_calc_model_t;
    typedef smart_ptr <well_rate_control_t, true>     sp_well_rate_control_t;

  public:

    MAKE_ME_HAPPY (well_rate_control_factory, objbase, "well_rate_control_factory");

    virtual sp_well_rate_control_t
    create_control (rate_control_type /*control_type*/, bool /*is_bhp*/, 
                    bool /*is_production*/, const sp_calc_model_t & /*calc_model*/)
    {
      // we have a problem with exporting pure abstract classes to puthon. 
      bs_throw_exception ("PURE_CALL"); 
    }

    virtual ~well_rate_control_factory ()
    {
    }
  };

} // namespace wells
} // namespace blue_sky


#endif  // #ifndef BS_WELLS_WELL_RATE_CONTROL_INTERFACE_H_

