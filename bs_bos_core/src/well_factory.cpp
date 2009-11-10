/**
 * \file well_factory.cpp
 * \brief impl of well and connection factories
 * \author Sergey Miryanov
 * \date 17.07.2008
 * */
#include "stdafx.h"

#include "calc_well.h"
#include "well_connection.h"
#include "reservoir.h"
#include "facility_manager.h"
#include "default_connection.h"
#include "default_well.h"

namespace blue_sky
  {

  template <typename strategy_t>
  well_factory<strategy_t>::well_factory (bs_type_ctor_param param /* = NULL */)
  {

  }
  template <typename strategy_t>
  well_factory<strategy_t>::well_factory (const well_factory<strategy_t> &f)
  : bs_refcounter (f), objbase (f)
  {

  }

  template <typename strategy_t>
  typename well_factory<strategy_t>::sp_well_t
  well_factory<strategy_t>::create_well (const std::string &group_name, const std::string &well_name) const
    {
      BS_ASSERT (well_name.size ());

      sp_well_t w = BS_KERNEL.create_object (wells::default_well <strategy_t>::bs_type (), true);
      BS_ASSERT (w) (well_name);

      w->set_name (well_name);

      //sp_well_t w = BS_KERNEL.create_object_v2 <wells::default_well <strategy_t> > (well_name, true);
      //BS_ERROR (w, "well_factory::create_well: Can't create well");

      return w;
    }

  template <typename strategy_t>
  typename well_factory<strategy_t>::sp_connection_t
  well_factory<strategy_t>::create_connection () const
    {
      sp_connection_t con = BS_KERNEL.create_object (wells::default_connection <strategy_t>::bs_type (), true);
      BS_ASSERT (con);

      return con;
    }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (well_factory, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (well_factory, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_factory <base_strategy_fi>), 1, (objbase), "well_factory_seq_fi", "Base class for well factory", "Base class for well factory", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_factory <base_strategy_di>), 1, (objbase), "well_factory_seq_di", "Base class for well factory", "Base class for well factory", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_factory <base_strategy_mixi>), 1, (objbase), "well_factory_seq_mixi", "Base class for well factory", "Base class for well factory", false);
  //////////////////////////////////////////////////////////////////////////

  bool
  well_factory_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, well_factory<base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, well_factory<base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, well_factory<base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);

    return res;
  }

}	// namespace blue_sky
