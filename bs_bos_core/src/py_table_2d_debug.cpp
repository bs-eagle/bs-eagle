/**
 *       \file  py_table_2d_debug.cpp
 *      \brief  Python wrappers for table_2d
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.05.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#include "stdafx.h"

#ifdef _DEBUG

#include "py_table_2d_debug.h"

namespace blue_sky
  {
  namespace python
    {

    using namespace boost::python;

    void export_py_data_row (const char *name)
    {
      class_ < py_data_row > (name)
      .def ("get", &py_data_row::get)
      ;
    }

    void export_py_data_group (const char *name)
    {
      class_ < py_data_group > (name)
      .def ("get_rows_count", &py_data_group::get_rows_count)
      .def ("get_columns_count", &py_data_group::get_columns_count)
      .def ("get_row", &py_data_group::get_row)
      ;
    }

    void export_py_table_2d (const char *name)
    {
      class_ < py_table_2d> (name)
      .def ("get_groups_count", &py_table_2d ::get_groups_count)
      .def ("get_data_group", &py_table_2d ::get_data_group)
      ;
    }

    void
    py_export_table_2d ()
    {

      export_py_data_row ("data_row");
      export_py_data_group ("data_group");
      export_py_table_2d ("table_2d");
    }


  }
}


#endif  // #ifdef _DEBUG
