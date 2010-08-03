/**
 * \file py_save_shared_vector.cpp
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

    struct py_save_shared_vector : public tools::save_shared_vector
      {

        py_save_shared_vector (const char *filename)
            : tools::save_shared_vector (filename)
        {

        }

        void save_list_int (shared_vector<int> &array)
        {
          save (array);
        }
        void save_list_float (shared_vector<float> &array)
        {
          save (array);
        }
        void save_list_double (shared_vector<double> &array)
        {
          save (array);
        }
      };

    void
    py_export_save_shared_vector ()
    {
      using namespace boost::python;

      class_ <py_save_shared_vector> ("save_shared_vector", init <const char *> ())
      .def ("save", &py_save_shared_vector::save_list_int)
      .def ("save", &py_save_shared_vector::save_list_float)
      .def ("save", &py_save_shared_vector::save_list_double)
      ;
    }

  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
