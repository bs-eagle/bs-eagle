/**
 *       \file  py_pvt.h
 *      \brief  python wrappers for pvt_
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.05.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef PY_BS_PVT_H_
#define PY_BS_PVT_H_

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky  {
namespace python    {

  void
  py_export_pvt ();

} // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
#endif  // #ifndef PY_BS_PVT_H_
