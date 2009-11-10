/**
 * \file py_naive_file_reader.cpp
 * \brief impl of
 * \author Sergey Miryanov
 * \date 18.06.2008
 *
 */
#include "bs_matrix_stdafx.h"

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

        void read_list_int (base_strategy_di::index_array_t &array)
        {
          read_list (array);
        }
        void read_list_float (base_strategy_fi::item_array_t &array)
        {
          read_list (array);
        }
        void read_list_double (base_strategy_di::item_array_t &array)
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
      ;
    }

  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
