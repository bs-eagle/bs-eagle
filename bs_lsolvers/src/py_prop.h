#ifndef PY_PROP_H_
#define PY_PROP_H_
/** 
 * @file py_prop.h
 * @brief 
 * @author 
 * @date 2009-11-25
 */
#include <string>
#include "prop_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

    // Sergey Miryanov at 07.04.2008
    // Refactored at 16.10.2009
    // for casting into python from child classes to parent class
    // we should export base class
    
  PY_EXPORTER (py_prop_exporter, default_exporter)
    .def ("add_property_f",                     &T::add_property_f, 
        args ("default_value", "name", "description"), "Add fp type property")
    .def ("add_property_i",                     &T::add_property_i, 
        args ("default_value", "name", "description"), "Add integer type property")
    .def ("add_property_s",                     &T::add_property_s, 
        args ("default_value", "name", "description"), "Add string type property")
    .def ("add_property_b",                     &T::add_property_b, 
        args ("default_value", "name", "description"), "Add boolean type property")
    .def ("clear",                              &T::clear, 
        args (""), "Clear all")
    .def ("get_index_f",                        &T::get_index_f, 
        args ("property_name"), "Return fp type property index by name")
    .def ("get_index_i",                        &T::get_index_i, 
        args ("property_name"), "Return int type property index by name")
    .def ("get_index_s",                        &T::get_index_s, 
        args ("property_name"), "Return string type property index by name")
    .def ("get_index_b",                        &T::get_index_b, 
        args ("property_name"), "Return boolean type property index by name")
    .def ("get_f",                              &T::get_f, 
        args ("property_index"), "Return fp type property value by index")
    .def ("get_i",                              &T::get_i, 
        args ("property_index"), "Return int type property value by index")
    .def ("get_s",                              &T::get_s, 
        args ("property_index"), "Return string type property value by index")
    .def ("get_b",                              &T::get_b, 
        args ("property_index"), "Return boolean type property value by index")
    .def ("set_f",                              &T::set_f,
        args ("property_index", "value"), "Set value by property index")
    .def ("set_i",                              &T::set_i,
        args ("property_index", "value"), "Set value by property index")
    .def ("set_s",                              &T::set_s,
        args ("property_index", "value"), "Set value by property index")
    .def ("set_b",                              &T::set_b,
        args ("property_index", "value"), "Set value by property index")
    .def ("check_f",                            &T::check_f,
        args ("property_index"), "Return false if property value set to default")
    .def ("check_i",                            &T::check_i,
        args ("property_index"), "Return false if property value set to default")
    .def ("check_s",                            &T::check_s,
        args ("property_index"), "Return false if property value set to default")
    .def ("check_b",                            &T::check_b,
        args ("property_index"), "Return false if property value set to default")
    .def ("reset_f",                            &T::reset_f,
        args ("property_index"), "Reset property value to default")
    .def ("reset_i",                            &T::reset_i,
        args ("property_index"), "Reset property value to default")
    .def ("reset_s",                            &T::reset_s,
        args ("property_index"), "Reset property value to default")
    .def ("reset_b",                            &T::reset_b,
        args ("property_index"), "Reset property value to default")
    .def ("reset_all",                          &T::reset_all,
        args (""), "Reset all properties to default")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  PY_EXPORTER (py_dummy_exporter, default_exporter)
                                                 
  PY_EXPORTER_END;                               
                                                 
  //! export matrices to python                  
  void py_export_prop ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_IFACE_H_
