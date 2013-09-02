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

  well_factory::well_factory (bs_type_ctor_param /*param  = NULL */)
  {

  }
  well_factory::well_factory (const well_factory &f)
  : bs_refcounter (f), objbase (f)
  {

  }

  well_factory::sp_well_t
  well_factory::create_well (const std::string & /*group_name*/, const std::string &well_name) const
    {
      BS_ASSERT (well_name.size ());

      sp_well_t w = BS_KERNEL.create_object (wells::default_well::bs_type (), true);
      BS_ASSERT (w) (well_name);

      w->set_name (well_name);

      //sp_well_t w = BS_KERNEL.create_object_v2 <wells::default_well > (well_name, true);
      //BS_ERROR (w, "well_factory::create_well: Can't create well");

      return w;
    }

  well_factory::sp_connection_t
  well_factory::create_connection () const
    {
      sp_connection_t con = BS_KERNEL.create_object (wells::default_connection::bs_type (), true);
      BS_ASSERT (con);

      return con;
    }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (well_factory);
  BLUE_SKY_TYPE_STD_COPY (well_factory);
  BLUE_SKY_TYPE_IMPL (well_factory, objbase, "well_factory", "Base class for well factory", "Base class for well factory");
  //////////////////////////////////////////////////////////////////////////

  bool
  well_factory_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, well_factory::bs_type ());
    BS_ASSERT (res);

    return res;
  }

}	// namespace blue_sky
