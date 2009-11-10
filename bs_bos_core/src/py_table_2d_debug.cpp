/**
 * \file py_table_2d_debug.cpp
 * \brief python wrappers for table_2d. debug only
 * \author Miryanov Sergey
 * \date 12.05.2008
 */
#include "stdafx.h"

#ifdef _DEBUG

#include "py_table_2d_debug.h"

namespace blue_sky
  {
  namespace python
    {

    using namespace boost::python;

    template <typename strategy_t>
    void export_py_data_row (const char *name)
    {
      class_ < py_data_row <strategy_t> > (name)
      .def ("get", &py_data_row <strategy_t>::get)
      ;
    }

    template <typename strategy_t>
    void export_py_data_group (const char *name)
    {
      class_ < py_data_group <strategy_t> > (name)
      .def ("get_rows_count", &py_data_group <strategy_t>::get_rows_count)
      .def ("get_columns_count", &py_data_group <strategy_t>::get_columns_count)
      .def ("get_row", &py_data_group <strategy_t>::get_row)
      ;
    }

    template <typename strategy_t>
    void export_py_table_2d (const char *name)
    {
      class_ < py_table_2d <strategy_t> > (name)
      .def ("get_groups_count", &py_table_2d <strategy_t>::get_groups_count)
      .def ("get_data_group", &py_table_2d <strategy_t>::get_data_group)
      ;
    }

    void
    py_export_table_2d ()
    {

      export_py_data_row <base_strategy_fi> ("data_row_fi");
      export_py_data_row <base_strategy_di> ("data_row_di");

      export_py_data_group <base_strategy_fi> ("data_group_fi");
      export_py_data_group <base_strategy_di> ("data_group_di");

      export_py_table_2d <base_strategy_fi> ("table_2d_fi");
      export_py_table_2d <base_strategy_di> ("table_2d_di");
    }


  }
}


#endif  // #ifdef _DEBUG
