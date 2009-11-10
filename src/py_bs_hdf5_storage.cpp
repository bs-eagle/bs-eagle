/**
 * \file py_bs_hdf5_storage.cpp
 * \brief Python wrappers for
 * \author
 * \date
 */
//#include "stdafx.h"
#include "bs_hdf5_storage.h"
#include "py_bs_hdf5_storage.h"

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky
  {
  namespace python
    {

    //! export wrappers
    void
    py_export_bs_hdf5_storage ()
    {

      using namespace boost::python;

      class_ <py_bs_hdf5_storage, bases <py_objbase> > ("hdf5_storage")
      .def ("open", &py_bs_hdf5_storage::open)
      .def ("close", &py_bs_hdf5_storage::close)
      .def ("writei", &py_bs_hdf5_storage::write<int>)
      .def ("writef", &py_bs_hdf5_storage::write<float>)
      .def ("writed", &py_bs_hdf5_storage::write<double>)
      ;
    }
  }
}

#endif // #ifdef BSPY_EXPORTING_PLUGIN

