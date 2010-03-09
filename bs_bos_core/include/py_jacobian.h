/**
 *       \file  py_jacobian.h
 *      \brief  Exports python wrapper for jacobian,
 *              see jacobian.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_PY_JACOBIAN_H_
#define BS_PY_JACOBIAN_H_

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky {
namespace python {

  /**
   * \brief  Exports wrapper to python
   * */
  void
  py_export_jacobian ();

} // namespace python
} // namespace blue_sky

#endif
#endif // #ifndef BS_PY_JACOBIAN_H_

