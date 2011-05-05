/**
 *       \file  main.cpp
 *      \brief  
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  05.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "std_well_stdafx.hpp"

#include "well_events.h"
#include "std_well_keywords.hpp"

using namespace blue_sky;
using namespace boost::python;

namespace 
{
  bool
  register_types (plugin_descriptor const &pd)
  {
    bool res = true;
    res &= BS_KERNEL.register_type (pd, WELSPECS_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WELLCON_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, COMPDAT_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WCONPROD_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WCONHIST_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WCONINJE_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WECON_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WECONINJ_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WEFAC_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WELTARG_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WPIMULT_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, COMPENSATION_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, PERMFRAC_event::bs_type ()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, std_well_keywords::bs_type ()); BS_ASSERT (res);

    return res;
  }
#ifdef BSPY_EXPORTING_PLUGIN
  void
  init_py_subsystem ()
  {
  }
#endif
}

namespace blue_sky 
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("std_well", "1.0.0", "STD Well", "STD Well", "std_well")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN 
{
  init_py_subsystem ();
}


#ifdef _DEBUG
BOOST_PYTHON_MODULE (std_well_d)
#else
BOOST_PYTHON_MODULE (std_well)
#endif
{
  init_py_subsystem ();
  bool res = register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    bs_throw_exception ("Can't register std_well types");
}

#endif
