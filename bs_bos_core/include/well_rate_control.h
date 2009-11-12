/**
 *       \file  well_rate_control.h
 *      \brief  Implementation of well control
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_H_
#define BS_WELLS_WELL_RATE_CONTROL_H_

#include "well_rate_control_interface.h"
#include "rate_control_type.h"

namespace blue_sky
  {
  namespace wells
    {

      template <typename strategy_t>
      class well_rate_control_interface;

    template <typename strategy_t>
    class BS_API_PLUGIN well_rate_control : public objbase
      {
      public:
        typedef strategy_t                                  strategy_type;
        typedef well_rate_control <strategy_t>              this_t;

        typedef typename strategy_t::index_t                index_t;

        typedef calc_model <strategy_t>                     calc_model_t;
        typedef jacobian_matrix <strategy_t>                jmatrix_t;
        typedef well <strategy_t>                           well_t;
        typedef well_controller <strategy_t>                well_controller_t;
        typedef well_rate_control_interface <strategy_t>    well_rate_control_impl_t;

        typedef smart_ptr <calc_model_t, true>              sp_calc_model_t;
        typedef smart_ptr <jmatrix_t, true>                 sp_jmatrix_t;
        typedef smart_ptr <well_t, true>                    sp_well_t;
        typedef smart_ptr <well_controller_t, true>         sp_well_controller_t;

        typedef smart_ptr <well_rate_control_impl_t, true>  sp_well_rate_control_impl_t;

      public:

        bool
        is_bhp () const
        {
          return is_bhp_;
        }
        bool
        is_production () const
        {
          return is_prod_;
        }
        bool
        is_rate () const
        {
          return !is_bhp_;
        }
        rate_control_type
        get_control_type () const
        {
          return control_type_;
        }

        void
        set_is_bhp (bool f)
        {
          is_bhp_ = f;
        }
        void 
        set_is_prod (bool f)
        {
          is_prod_ = f;
        }
        void 
        set_control_type (rate_control_type control_type)
        {
          control_type_ = control_type;
        }
        void
        set_impl (const sp_well_rate_control_impl_t &impl)
        {
          impl_ = impl;
        }
        void 
        compute_rate (const sp_calc_model_t &calc_model, sp_jmatrix_t &jmatrix, sp_well_t &well, const sp_well_controller_t &well_controller) const
        {
          impl_->compute_rate (calc_model, jmatrix, well, well_controller);
        }
        void 
        compute_derivs (const sp_calc_model_t &calc_model, sp_jmatrix_t &jmatrix, sp_well_t &well, const sp_well_controller_t &well_controller) const
        {
          impl_->compute_derivs (calc_model, jmatrix, well, well_controller);
        }

        BLUE_SKY_TYPE_DECL_T (well_rate_control);

      private:
        bool                                            is_bhp_;
        bool                                            is_prod_;
        auto_value <rate_control_type, null_control>    control_type_;
        sp_well_rate_control_impl_t                     impl_;
      };

  } // namespace wells
} // namespace blue_sky

#endif // #ifndef BS_WELLS_WELL_RATE_CONTROL_H_

