/** 
 * @file py_table.h
 * @brief Python #table and #table_iface interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */
#ifndef PY_TABLE_PZHJ8KVY

#define PY_TABLE_PZHJ8KVY

#include <string>
#include "table_iface.h"
#include "table.h"

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
#if 0 
      struct table_pickle_suite : boost::python::pickle_suite
        {
          static
          boost::python::tuple
          getinitargs(const table &w)
          {
            return boost::python::make_tuple();
          }

#if 0         
          static boost::python::tuple getstate (const table &w)
          {
            std::string &str = w.to_str ();
            return boost::python::make_tuple (str);
          }
          static void setstate(table &w, boost::python::tuple state)
          {
            using namespace boost::python;
            if (len(state) != 1)
            {
              PyErr_SetObject(PyExc_ValueError,
                ("expected 1-item tuple in call to __setstate__; got %s"
                % state).ptr()
                );
              throw_error_already_set();
            }

            std::string str = extract<std::string>(state[0]);
            w.from_str (str);

          }
#else
          static boost::python::tuple getstate (boost::python::object w_obj)
          {
            using namespace boost::python;
            table const& w = extract<table const&>(w_obj)();
            dict d = extract<dict>(w_obj.attr("__dict__"))();
            //std::string prj_file = extract<std::string>(w_obj.attr ("project_filename") );
            /*
            list iterkeys = (list) d.keys();
            for (int i = 0; i < len(iterkeys); i++)
              {
                std::string key = extract<std::string>(iterkeys[i]);
                if (key == "project_filename")
                  {
                    std::string value = extract<std::string>(d[iterkeys[i]]);
                  }
              }
            */
            std::string &str = w.to_str ();
            return make_tuple (w_obj.attr("__dict__"), str);
          }

          static
          void
          setstate(boost::python::object w_obj, boost::python::tuple state)
          {
            using namespace boost::python;
            table& w = extract<table&>(w_obj)();

            if (len(state) != 2)
            {
              PyErr_SetObject(PyExc_ValueError,
                ("expected 2-item tuple in call to __setstate__; got %s"
                % state).ptr()
                );
              throw_error_already_set();
            }

            // restore the object's __dict__
            dict d = extract<dict>(w_obj.attr("__dict__"))();
            d.update(state[0]);

            // restore the internal state of the C++ object
            std::string str = extract<std::string>(state[1]);
            w.from_str (str);
          }

          static bool getstate_manages_dict() { return true; }
#endif

        };
#endif 
      //using namespace ::boost::python;
      //class_<table_iface>("table_iface")
        // ...
        //.def_pickle(table_pickle_suite ());
        // ...

  PY_EXPORTER (py_table_exporter, default_exporter)
    .def ("init",                               &T::init, 
        args ("n_rows", "n_cols"), "Initialize table n_rows * n_cols")
    .def ("clear",                               &T::clear, 
        args (""), "Clear all")
    .def ("copy",                               &T::copy, 
        args ("table"), "Copy from another table")
    .def ("set_col_name",                       &T::set_col_name, 
        args ("col_idx", "name"), "Set name for the col_idx column")
    .def ("get_col_name",                       &T::get_col_name, 
        args ("col_idx"), "Return name of the col_idx column")
    .def ("get_n_rows",                         &T::get_n_rows, 
        args (""), "Return number of rows in table")
    .def ("get_n_cols",                         &T::get_n_cols, 
        args (""), "Return number of columns in table")
    .def ("set_col_values",                     &T::set_col_values, 
        args ("col_idx", "column"), "Set values of the column col_idx")
    .def ("get_col_values",                     &T::get_col_values, 
        args ("col_idx"), "Return values of the column col_idx")
    .def ("add_row",                            &T::add_row, 
        args ("row_idx"), "Insert row at the position row")
    .def ("remove_row",                            &T::remove_row, 
        args ("row_idx"), "Remove row at the position row")
    .def ("get_value",                            &T::get_value, 
        args ("row_idx", "column_idx"), "Return value at given row and column")
    .def ("set_value",                            &T::set_value, 
        args ("row_idx", "column_idx", "value"), "Set value at given row and column")
    .def ("check_serial",                         &T::check_serial, 
        args (""), "should return copy of the table")
    .def ("to_str",                               &T::to_str, 
        args (""), "serialize class content to string")
    .def ("from_str",                             &T::from_str, 
        args ("string"), "restore class content from string")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  //PY_EXPORTER (py_dummy_exporter, default_exporter)
                                                 
  //PY_EXPORTER_END;                               
                                                 
  //! export matrices to python                  
  void py_export_table ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN


#endif /* end of include guard: PY_TABLE_PZHJ8KVY */
