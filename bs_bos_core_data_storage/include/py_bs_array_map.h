#ifndef __PY_POOL_H
#define __PY_POOL_H
/**
 * \file   py_bs_pool.h
 * \brief  Python wrapper for pool template class
 * \author Mark Khait
 * \date 2010-08-26
 */

#include "bs_array_map.h"
#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {
      PY_EXPORTER (bs_array_map_exporter, default_exporter)
        .def("init",
            &T::init,
            args ("nx", "ny", "nz"), "Initialize bs_array_map by model dimensions")
        .def("add_item",
            &T::add_item,
            args ("key", "array", "dimens","def_val"), "Add array item")
        .def("create_item",
            &T::py_create_item,
            args ("key", "dimens","def_val"), "Create array item")
        .def("list_items",
            &T::py_list_items,
            args (""), "List all items")
      
      //.add_property ("title",               &T::title)
      PY_EXPORTER_END;

      
      void py_export_array_maps();
    } // namespace python
  } // namespace blue_sky

#endif // BSPY_EXPORTING_PLUGIN
#endif // __PY_POOL_H

