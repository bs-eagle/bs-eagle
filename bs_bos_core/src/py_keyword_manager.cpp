/**
 * \file py_keyword_manager.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 28.10.2009
 * */
#include "stdafx.h"
#include "py_keyword_manager.h"
#include "keyword_manager.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "read_class.h"
#include "data_class.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

namespace blue_sky {
namespace python {

  template <typename strategy_t>
  struct py_keyword_handler_iface_base : keyword_handler_iface <strategy_t>
  {
    typedef keyword_params <strategy_t> keyword_params_t;

    void 
    handler (const std::string &, keyword_params_t &)
    {
      bs_throw_exception ("PURE CALL FROM PYTHON");
    }
  };
  template <typename strategy_t>
  struct py_keyword_handler_iface : py_keyword_handler_iface_base <strategy_t>, boost::python::wrapper <py_keyword_handler_iface_base <strategy_t> >
  {
    typedef keyword_params <strategy_t> keyword_params_t;

    py_keyword_handler_iface ()
    {
    }
    py_keyword_handler_iface (const py_keyword_handler_iface_base <strategy_t> &)
    {
    }

    WRAP_PURE_METHOD (handler, void, 2, (const std::string &, keyword_params_t &));
  };

  template <typename T>
  void
  register_keyword (T *t, const std::string &keyword, const typename T::shared_handler_t &handler)
  {
    t->register_keyword (keyword, handler);
  }

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

  template <typename T>
  smart_ptr <idata, true>
  get_data (T *t)
  {
    typedef smart_ptr <idata, true> sp_idata_t;
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

  PY_EXPORTER (keyword_manager_exporter, default_exporter)
    .def ("register_keyword", register_keyword <T>)
    .def ("is_keyword_supported", &T::is_keyword_supported)
    .def ("is_keyword_activated", &T::is_keyword_activated)
  PY_EXPORTER_END;

  PY_EXPORTER (keyword_params_exporter, empty_exporter)
    .add_property ("reader", make_function (get_reader <T>))
    .add_property ("data",   make_function (get_data <T>))
    .add_property ("mesh",   make_function (get_mesh <T>))
    .add_property ("keyword_manager", make_function (get_keyword_manager <T>))
  PY_EXPORTER_END;

  PY_EXPORTER (keyword_handler_iface_exporter, empty_exporter)
    .def ("handler", &T::handler)
  PY_EXPORTER_END;

  void
  export_keyword_manager ()
  {
    using namespace boost::python;

    strategy_exporter::export_base_ext <keyword_params, keyword_params_exporter, class_type::concrete_class> ("keyword_params");
    strategy_exporter::export_base_ext <keyword_handler_iface, empty_exporter, class_type::abstract_class> ("keyword_handler_iface");
    strategy_exporter::export_class_ext <py_keyword_handler_iface, keyword_handler_iface, keyword_handler_iface_exporter, class_type::concrete_class> ("py_keyword_handler_iface");
    strategy_exporter::export_base <keyword_manager, keyword_manager_exporter> ("keyword_manager");
  }

} // namespace python
} // namespace blue_sky

