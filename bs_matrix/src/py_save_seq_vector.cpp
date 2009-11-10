/**
 * \file py_save_seq_vector.cpp
 * \brief impl of
 * \author Sergey Miryanov
 * \date 20.06.2008
 */
#include "bs_matrix_stdafx.h"

#include "py_save_seq_vector.h"
#include "save_seq_vector.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

    struct py_save_seq_vector : public tools::save_seq_vector
      {

        py_save_seq_vector (const char *filename)
            : tools::save_seq_vector (filename)
        {

        }

        void save_list_int (base_strategy_di::index_array_t &array)
        {
          save (array);
        }
        void save_list_float (base_strategy_fi::item_array_t &array)
        {
          save (array);
        }
        void save_list_double (base_strategy_di::item_array_t &array)
        {
          save (array);
        }
      };

    void
    py_export_save_seq_vector ()
    {
      using namespace boost::python;

      class_ <py_save_seq_vector> ("save_seq_vector", init <const char *> ())
      .def ("save", &py_save_seq_vector::save_list_int)
      .def ("save", &py_save_seq_vector::save_list_float)
      .def ("save", &py_save_seq_vector::save_list_double)
      ;
    }

  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
