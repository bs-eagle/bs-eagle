/**
 * \file py_pvt.cpp
 * \brief python wrappers
 * \author Miryanov Sergey
 * \date 08.05.2008
 */
#include "bs_pvt_stdafx.h"

#include "py_pvt.h"

#ifdef BSPY_EXPORTING_PLUGIN

#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "pvt_oil.h"
#include "pvt_water.h"

namespace blue_sky
  {
  namespace python
    {

    template <class pvt_t> void
    py_export (const char *name)
    {
      using namespace boost::python;

      class_<pvt_t, bases <typename pvt_t::base_t> > (name)
      .def ("build", &pvt_t::build)
      .def ("get_p_step", &pvt_t::get_p_step)
      .def ("get_surface_density", &pvt_t::get_surface_density)
      .def ("set_surface_density", &pvt_t::set_surface_density)
      ;
    }

    template <typename pvt_t>
    void
    py_pvt_export (const char * name)
    {
      using namespace boost::python;

      class_ <pvt_t, bases <typename pvt_t::base_t> > (name, init <const sp_obj &> ())
        .def ("calc", &pvt_t::calc)
        ;
    }

#define PY_EXPORT_PVT_WATER \
    py_export <py_pvt_wrapper <pvt_water <base_strategy_fi> > > ("pvt_water_fi_base");  \
    py_export <py_pvt_wrapper <pvt_water <base_strategy_di> > > ("pvt_water_di_base");  \
    py_pvt_export <py_pvt_water <pvt_water <base_strategy_fi> > > ("pvt_water_fi");   \
    py_pvt_export <py_pvt_water <pvt_water <base_strategy_di> > > ("pvt_water_di");

#define PY_EXPORT_PVT_GAS \
    py_export <py_pvt_wrapper <pvt_gas <base_strategy_fi> > > ("pvt_gas_fi_base");      \
    py_export <py_pvt_wrapper <pvt_gas <base_strategy_di> > > ("pvt_gas_di_base");      \
    py_pvt_export <py_pvt_gas <pvt_gas <base_strategy_fi> > > ("pvt_gas_fi");           \
    py_pvt_export <py_pvt_gas <pvt_gas <base_strategy_di> > > ("pvt_gas_di");

#define PY_EXPORT_PVT_OIL \
    py_export <py_pvt_wrapper <pvt_oil <base_strategy_fi> > > ("pvt_oil_fi_base");      \
    py_export <py_pvt_wrapper <pvt_oil <base_strategy_di> > > ("pvt_oil_di_base");      \
    py_pvt_export <py_pvt_oil <pvt_oil <base_strategy_fi> > > ("pvt_oil_fi");           \
    py_pvt_export <py_pvt_oil <pvt_oil <base_strategy_di> > > ("pvt_oil_di");

    void
    py_export_pvt ()
    {
      using namespace boost::python;

      class_ <py_pvt_base, bases <py_objbase> > ("pvt_base", init <sp_obj> ());

      PY_EXPORT_PVT_WATER;
      PY_EXPORT_PVT_GAS;
      PY_EXPORT_PVT_OIL;

      py_export<py_pvt_wrapper< pvt_dead_oil<base_strategy_fi> > >  ("pvt_dead_oil_fi");
      py_export<py_pvt_wrapper< pvt_dead_oil<base_strategy_di> > >  ("pvt_dead_oil_di");
    }

  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN

