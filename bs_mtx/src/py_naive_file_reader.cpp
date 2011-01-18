/**
 * \file py_naive_file_reader.cpp
 * \brief impl of
 * \author Sergey Miryanov
 * \date 18.06.2008
 *
 */
#include "bs_mtx_stdafx.h"

#include "py_naive_file_reader.h"
#include "naive_file_reader.h"

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky
  {
  namespace python
    {

    struct py_naive_file_reader : public naive_file_reader
      {

        py_naive_file_reader (const char *filename)
            : naive_file_reader (filename)
        {

        }
        void read_item_int (int &item)
        {
          read_item (item);
        }
        
        void read_list_int (shared_vector<int> &array)
        {
          read_list (array);
        }
        void read_list_float (shared_vector<float> &array)
        {
          read_list (array);
        }
        void read_list_double (shared_vector<double> &array)
        {
          read_list (array);
        }
        
      };

    void
    py_export_naive_file_reader ()
    {
      using namespace boost::python;

      class_ <py_naive_file_reader> ("naive_file_reader", init <const char *> ())
      .def ("read_list", &py_naive_file_reader::read_list_int)
      .def ("read_list", &py_naive_file_reader::read_list_float)
      .def ("read_list", &py_naive_file_reader::read_list_double)
      .def ("read_item", &py_naive_file_reader::read_item_int)
      .def ("locate_section", &py_naive_file_reader::locate_section, return_value_policy <reference_existing_object> ())
      ;
    }

  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
