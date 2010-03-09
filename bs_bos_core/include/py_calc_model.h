/**
 *       \file  py_calc_model.h
 *      \brief  Python wrappers for calc_model, calc_model_data
 *     \author  Sergey Miryanov
 *       \date  02.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef PY_CALC_MODEL_H
#define PY_CALC_MODEL_H

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky {
namespace python  {

  /**
   * \brief  Exports wrappers to python
   * */
  void 
  py_export_calc_model ();

} // namespace python
} // namespace blue_sky

#endif
#endif // PY_CALC_MODEL_H
