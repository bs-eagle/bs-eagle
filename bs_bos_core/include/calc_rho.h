/**
 * \file calc_rho.h
 * \brief class for rho calculation
 * \author Sergey Miryanov
 * \date 26.09.2008
 * */
#ifndef BS_CALC_RHO_H_
#define BS_CALC_RHO_H_

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  class well;

  template <typename strategy_t>
  class calc_rho_base : public objbase
    {
    public:
      typedef typename strategy_t::index_t          index_t;
      typedef typename strategy_t::item_t           item_t;
      typedef typename strategy_t::item_array_t     item_array_t;

      typedef calc_model <strategy_t>               calc_model_t;
      typedef well <strategy_t>                     well_t;
      typedef rs_mesh_iface <strategy_t>                  mesh_iface_t;

      typedef smart_ptr <calc_model_t, true>        sp_calc_model_t;
      typedef smart_ptr <well_t, true>              sp_well_t;
      typedef smart_ptr <mesh_iface_t, true>              sp_mesh_iface_t;

    public:
      virtual ~calc_rho_base () {}

      virtual void calculate (const sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const = 0;
    };

  template <typename strategy_t>
  class calc_total_average_rho : public calc_rho_base <strategy_t>
    {
    public:
      typedef typename strategy_t::index_t          index_t;
      typedef typename strategy_t::item_t           item_t;
      typedef typename strategy_t::item_array_t     item_array_t;

      typedef calc_model <strategy_t>               calc_model_t;
      typedef well <strategy_t>                     well_t;
      typedef rs_mesh_iface <strategy_t>                  mesh_iface_t;

      typedef smart_ptr <calc_model_t, true>        sp_calc_model_t;
      typedef smart_ptr <well_t, true>              sp_well_t;
      typedef smart_ptr <mesh_iface_t, true>              sp_mesh_iface_t;

    public:

      void calculate (const sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const;

      BLUE_SKY_TYPE_DECL_T (calc_total_average_rho <strategy_t>);
    };

  bool
  calc_rho_register_types (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky

#endif  // #ifndef BS_CALC_RHO_H_

