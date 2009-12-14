/**
 *       \file  py_scal_wrapper.h
 *      \brief  Exports scal_3p (+) to python
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  14.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef PY_BS_SCAL_WRAPPER_H_
#define PY_BS_SCAL_WRAPPER_H_

#include "scal_3p.h"
#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"
#include "jfunction.h"

namespace blue_sky  {
namespace python    {

  void py_export_scal ();

} // namespace python
} // namespace blue_sky


#endif // #ifndef PY_BS_SCAL_WRAPPER_H_
