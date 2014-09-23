/**
 *       \file  py_data_diments.cpp
 *      \brief  Implements export of data_dimens to python
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.11.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#ifdef BSPY_EXPORTING_PLUGIN
#include "data_dimens.h"
#include <boost/python.hpp>

namespace blue_sky {
namespace python {

  void
  export_data_dimens ()
  {
    using namespace boost::python;

    class_ <data_dimens> ("dimens", init < t_long, t_long, t_long > ())
      .add_property ("nx", &data_dimens::nx)
      .add_property ("ny", &data_dimens::ny)
      .add_property ("nz", &data_dimens::nz)
      ;
  }

} // namespace python
} // namespace blue_sky
#endif
