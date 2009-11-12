/**
 *       \file  simulator_events.h
 *      \brief  Tools to declare events which will be
 *              called in main calculation loop
 *    \details  DECLARE_EVENT_LIST declares by sequence of declaration 
 *              a set of functions with prefix on_, 
 *              see reservoir_simulator.h for usage example
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_SIMULATOR_EVENTS_H_
#define BS_SIMULATOR_EVENTS_H_

#include "pp_param_list.h"


#define PUSH_LIST(seq_, size_) BOOST_PP_CAT(PUSH_LIST_, size_) seq_
#define PUSH_LIST_0
#define PUSH_LIST_1(x) params->push_back(boost::lexical_cast<std::string> (A1));
#define PUSH_LIST_2(x) params->push_back(boost::lexical_cast<std::string> (A2)); PUSH_LIST_1
#define PUSH_LIST_3(x) params->push_back(boost::lexical_cast<std::string> (A3)); PUSH_LIST_2
#define PUSH_LIST_4(x) params->push_back(boost::lexical_cast<std::string> (A4)); PUSH_LIST_3

#define DECL(type_, name_, t_, ts_)                                                                   \
  void BOOST_PP_CAT(on_, name_) (PARAM_LIST_FOR_EACH(t_, ts_))                                        \
  {                                                                                                   \
    smart_ptr <type_> params = BS_KERNEL.create_object (type_::bs_type ());                           \
    PUSH_LIST (BOOST_PP_TUPLE_TO_SEQ(ts_, t_), ts_)                                                   \
    this->fire_signal (BOOST_PP_CAT(name_, _signal), params);                                         \
  }

#define UDECL_I(r, data, i, elem)                                                                     \
  DECL(data, BOOST_PP_TUPLE_ELEM(3,0,elem), BOOST_PP_TUPLE_ELEM(3,1,elem), BOOST_PP_TUPLE_ELEM(3,2,elem))

#define DECLARE_SIGNAL_ENUM_I(r, data, i, elem)                                                       \
  BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(3,0,elem), _signal),

#define DECLARE_SIGNALS_ENUM(ts_, t_)                                                                 \
  BLUE_SKY_SIGNALS_DECL_BEGIN(bs_node)                                                                \
  BOOST_PP_SEQ_FOR_EACH_I(DECLARE_SIGNAL_ENUM_I, _, BOOST_PP_TUPLE_TO_SEQ(ts_, t_))                   \
  BLUE_SKY_SIGNALS_DECL_END

#ifdef BSPY_EXPORTING_PLUGIN
#define EXPORT_SIGNAL_ENUM_I(r, data, i, elem)                                                        \
  .value (BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(3,0,elem)), data::BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(3,0,elem), _signal))

#define EXPORT_SIGNALS_ENUM(name_, ts_, t_)                                                           \
  static void py_export_signals_enum ()                                                               \
  {                                                                                                   \
    using namespace boost::python;                                                                    \
    enum_ <signal_codes> (BOOST_PP_STRINGIZE(BOOST_PP_CAT(name_, _signals)))                          \
      BOOST_PP_SEQ_FOR_EACH_I(EXPORT_SIGNAL_ENUM_I, name_, BOOST_PP_TUPLE_TO_SEQ(ts_, t_))            \
      .export_values ()                                                                               \
      ;                                                                                               \
  }
#endif

#ifdef BSPY_EXPORTING_PLUGIN
#define DECLARE_EVENT_LIST(name_, params_, ts_, t_)                                                   \
  public:                                                                                             \
  DECLARE_SIGNALS_ENUM(ts_, t_)                                                                       \
  BOOST_PP_SEQ_FOR_EACH_I(UDECL_I, params_, BOOST_PP_TUPLE_TO_SEQ(ts_, t_))                           \
  EXPORT_SIGNALS_ENUM(name_, ts_, t_)
#else
#define DECLARE_EVENT_LIST(name_, params_, ts_, t_)                                                   \
  public:                                                                                             \
  DECLARE_SIGNALS_ENUM(ts_, t_)                                                                       \
  BOOST_PP_SEQ_FOR_EACH_I(UDECL_I, params_, BOOST_PP_TUPLE_TO_SEQ(ts_, t_))
#endif



#endif  // #ifndef BS_SIMULATOR_EVENTS_H_
