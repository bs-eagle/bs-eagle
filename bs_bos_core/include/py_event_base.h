/**
 *       \file  py_event_base.h
 *      \brief  Python wrappers for event_base, python events iterator,
 *              for event_base see event_base.h
 *     \author  Nikonov Max
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */
#ifndef PY_EVENT_BASE_H
#define PY_EVENT_BASE_H

#include "event_base.h"
#include "event_manager.h"

#include "calc_model.h"
#include "reservoir.h"

#include "export_python_wrapper.h"

#include "py_reservoir.h"
#include "py_calc_model.h"

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"

namespace blue_sky {
namespace python {

    STRATEGY_CLASS_WRAPPER (event_base, py_event_base)
    {
    public:
      typedef rs_mesh_iface <strategy_t>        mesh_t;
      typedef reservoir <strategy_t>            reservoir_t;
      typedef smart_ptr <mesh_t, true>          sp_mesh_iface_t;
      typedef smart_ptr <reservoir_t, true>     sp_reservoir_t;
      typedef event_base <strategy_t>           base_t;

    public:
      MAKE_ME_HAPPY (py_event_base, py_event_base_base <strategy_t>, "py_event_base");
      //WRAPPER_METHOD (apply, void, 3, (const sp_reservoir_t &, const py_mesh_iface_t &, const py_calc_model_t &));
    };

    /**
     * \brief  Exports wrappers to python
     * */
    void 
    py_export_events ();

} // namespace python
} // namespace blue_sky

#endif // PY_EVENT_BASE_H
