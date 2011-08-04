/**
 *       \file  py_keyword_manager.cpp
 *      \brief  Exports python wrappers for keyword_manager,
 *              see keyword_manager.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  28.10.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "bs_bos_core_data_storage_stdafx.h"
#include "py_keyword_manager.h"
#include "keyword_manager.h"
#include "read_class.h"
#include "data_class.h"
#include "rs_mesh_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()

#include BS_STOP_PLUGIN_IMPORT ()

#ifdef BSPY_EXPORTING_PLUGIN
#include "export_python_wrapper.h"

namespace blue_sky {
namespace python {

  
  struct py_keyword_handler_iface_base : keyword_handler_iface 
  {
    typedef keyword_params  keyword_params_t;

    void 
    handler (const std::string &, keyword_params_t &)
    {
      bs_throw_exception ("PURE CALL FROM PYTHON");
    }
  };
  
  struct py_keyword_handler_iface : py_keyword_handler_iface_base , boost::python::wrapper <py_keyword_handler_iface_base  >
  {
    typedef keyword_params  keyword_params_t;

    py_keyword_handler_iface ()
    {
    }
    py_keyword_handler_iface (const py_keyword_handler_iface_base  &)
    {
    }

    WRAP_PURE_METHOD (handler, void, 2, (const std::string &, keyword_params_t &));
  };

  //template <typename T>
  //void
  //register_keyword (T *t, const std::string &keyword, const typename T::shared_handler_t &handler, bool replace_existing)
  //{
  //  t->register_keyword (keyword, handler, replace_existing);
  //}

  template <typename T>
  smart_ptr <FRead, true>
  get_reader (T *t)
  {
    typedef smart_ptr <FRead, true> sp_reader_t;
    sp_reader_t reader (t->reader, bs_dynamic_cast ());

    if (!reader)
      {
        bs_throw_exception (boost::format ("Can't cast params->reader to sp_reader_t (%s)") % t->reader->bs_resolve_type ().stype_);
      }

    return reader;
  }

  /*
  template <typename T>
  smart_ptr <idata<typename T::strategy_type>, true>
  get_data (T *t)
  {
    typedef smart_ptr <idata<typename T::strategy_type>, true> sp_idata_t;
    sp_idata_t data (t->data, bs_dynamic_cast ());

    if (!data)
      {
        bs_throw_exception (boost::format ("Can't cast params->data to sp_idata_t (%s)") % t->data->bs_resolve_type ().stype_);
      }

    return data;
  }

  template <typename T>
  smart_ptr <rs_mesh_iface <typename T::strategy_type>, true>
  get_mesh (T *t)
  {
    typedef smart_ptr <rs_mesh_iface <typename T::strategy_type>, true> sp_mesh_t;
    sp_mesh_t mesh (t->mesh, bs_dynamic_cast ());

    if (!mesh)
      {
        bs_throw_exception (boost::format ("Can't cast params->mesh to sp_mesh_t (%s)") % t->mesh->bs_resolve_type ().stype_);
      }

    return mesh;
  }
  
  template <typename T>
  smart_ptr <keyword_manager <typename T::strategy_type>, true>
  get_keyword_manager (T *t)
  {
    typedef smart_ptr <keyword_manager <typename T::strategy_type>, true> sp_keyword_manager_t;
    sp_keyword_manager_t km (t->km, bs_dynamic_cast ());

    if (!km)
      {
        bs_throw_exception (boost::format ("Can't cast params->km to sp_keyword_manager_t (%s)") % t->km->bs_resolve_type ().stype_);
      }

    return km;
  }
  */

  PY_EXPORTER (keyword_manager_exporter, default_exporter)
    //.def ("register_keyword", register_keyword <T>)
    .def ("register_i_pool_keyword", &T::py_register_i_pool_keyword)
    .def ("register_fp_pool_keyword", &T::py_register_fp_pool_keyword)
    .def ("register_keywords", &T::register_plugin_keywords)
    .def ("register_plugin_keywords", &T::register_plugin_keywords)
    .def ("list_active_keywords", &T::py_list_active_keywords)
    .def ("list_supported_keywords", &T::py_list_supported_keywords)
    .def ("is_keyword_supported", &T::is_keyword_supported)
    .def ("is_keyword_activated", &T::is_keyword_activated)
    .def ("init", &T::init)
  PY_EXPORTER_END;

  PY_EXPORTER (keyword_params_exporter, empty_exporter)
    .add_property ("hdm", make_function (get_reader <T>))
    //.add_property ("data",   make_function (get_data <T>))
    //.add_property ("mesh",   make_function (get_mesh <T>))
    //.add_property ("keyword_manager", make_function (get_keyword_manager <T>))
  PY_EXPORTER_END;

  PY_EXPORTER (keyword_handler_iface_exporter, empty_exporter)
    .def ("handler", &T::handler)
  PY_EXPORTER_END;

  void
  export_keyword_manager ()
  {
    using namespace boost::python;

    base_exporter<keyword_manager_iface, empty_exporter>::export_class ("keyword_manager_iface");
    //class_exporter<keyword_manager, keyword_manager_iface, keyword_manager_exporter>::export_class ("keyword_manager");


    //strategy_exporter::export_base_ext <keyword_params, keyword_params_exporter, class_type::concrete_class> ("keyword_params");
	
    /*
    strategy_exporter::export_base_ext <keyword_handler_iface, empty_exporter, class_type::abstract_class> ("keyword_handler_iface");
    strategy_exporter::export_class_ext <py_keyword_handler_iface, keyword_handler_iface, keyword_handler_iface_exporter, class_type::concrete_class> ("py_keyword_handler_iface");
    strategy_exporter::export_base <keyword_manager, keyword_manager_exporter> ("keyword_manager");
    */
  }

} // namespace python
} // namespace blue_sky
#endif
