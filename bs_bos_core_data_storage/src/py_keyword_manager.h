/**
 *       \file  py_keyword_manager.h
 *      \brief  Exports python wrappers for keyword_manager,
 *              see keyword_manager.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  28.10.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_BOS_CORE_PY_KEYWORD_MANAGER_H_
#define BS_BOS_CORE_PY_KEYWORD_MANAGER_H_

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky {
namespace python {

  /**
   * \brief  Exports wrappers to python
   * */
  void
  export_keyword_manager ();

} // namespace python
} // namespace blue_sky

#endif
#endif // #ifndef BS_BOS_CORE_PY_KEYWORD_MANAGER_H_