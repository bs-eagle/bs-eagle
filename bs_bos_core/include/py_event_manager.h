/**
 *       \file  py_event_manager.h
 *      \brief  Python wrapper for event_manager, 
 *              see event_manager.h
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_PY_EVENT_MANAGER_H
#define BS_PY_EVENT_MANAGER_H

#ifdef BSPY_EXPORTING_PLUGIN

#include "py_event_base.h"
#include "event_manager.h"

namespace blue_sky {
namespace python {

  template <typename strategy_t>
  class event_creater_base
  {
  public:
    typedef py_event_base <strategy_t>  py_event_base_t;

    event_creater_base (boost::python::object &pobj);

    py_event_base_t *create ();

  protected:
    boost::python::object &obj;
  };

  /**
   * \brief Python wrapper to export event_manager to python
   * */
  template <typename strategy_t>
  class py_event_manager : public py_objbase
  {
  public:
    typedef py_event_manager                 this_t;
    typedef event_manager <strategy_t>       wrapped_t;
    typedef smart_ptr<wrapped_t, true>       sp_em_t;
    typedef py_event_base <strategy_t>       py_event_base_t;
    typedef event_creater_base <strategy_t>  event_creater_base_t;

    typedef event_list_iterator <strategy_t> event_list_iterator_t;

  public:
    py_event_manager();
    py_event_manager(const sp_em_t &src);
    py_event_manager(const this_t &src);

    ~py_event_manager ();

    py_event_base_t *create_event(const std::string &date, const std::string &event_name, const std::string &event_params);

    void factory_register(const std::string &name, const event_creater_base_t &creater);

    event_list_iterator_t el_begin ();
    //event_list_iterator_t el_end ();
  };


  /**
   * \brief  Exports wrappers to python
   * */
  void
  py_export_event_manager ();

} // namespace python
} // namespace blue_sky

#endif //#ifdef BSPY_EXPORTING_PLUGIN
#endif //#ifndef BS_PY_EVENT_MANAGER_H
