/**
 *       \file  well_event_params_decl.h
 *      \brief  Parameters declaration for well events
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  15.07.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_BOS_CORE_WELL_EVENT_PARAMS_DECL_H_
#define BS_BOS_CORE_WELL_EVENT_PARAMS_DECL_H_

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/stringize.hpp>

#include "named_pbase_type_helper.h"

#define PARAMS_DECL_BLUE_SKY(type_t, base_t, type_name)             \
  BLUE_SKY_TYPE_STD_CREATE_T_MEM (type_t)                           \
  BLUE_SKY_TYPE_STD_COPY_T_MEM (type_t)                             \
  BS_LOCK_THIS_DECL(type_t);                                        \
                                                                    \
  friend class type_descriptor;                                     \
                                                                    \
  static const type_descriptor &                                    \
  td_maker (const std::string &stype_postfix)                       \
  {                                                                 \
    static blue_sky::type_descriptor td(Loki::Type2Type<type_t> ()  \
      , Loki::Type2Type <base_t> ()                                 \
      , Loki::Int2Type <false> ()                                   \
      , stype_postfix                                               \
      , ""                                                          \
      , "");                                                        \
                                                                    \
    return td;                                                      \
  }                                                                 \
                                                                    \
  static blue_sky::type_descriptor bs_type()                        \
  {                                                                 \
    return td_maker (std::string (type_name) + "_" + BOOST_CURRENT_FUNCTION); \
  }                                                                 \
  virtual blue_sky::type_descriptor bs_resolve_type() const         \
  {                                                                 \
    return bs_type ();                                              \
  }


#define PARAMS_DECLARE_CLASS_I(r, data, i, elem)                              \
  BOOST_PP_TUPLE_ELEM (3, 0, elem),

#define PARAMS_DECLARE_CLASS_GETTERS_I(r, data, i, elem)                      \
  tools::named_pbase_value_type_helper <BOOST_PP_TUPLE_ELEM (3, 2, elem)>::type \
  BOOST_PP_CAT (get_, BOOST_PP_TUPLE_ELEM (3, 0, elem)) () const              \
  {                                                                           \
    idx_type index = BOOST_PP_TUPLE_ELEM (3, 0, elem);                        \
    return tools::named_pbase_value_type_helper <BOOST_PP_TUPLE_ELEM (3, 2, elem)>::get (this, index); \
  }                                                                           \
  tools::named_pbase_value_type_helper <BOOST_PP_TUPLE_ELEM (3, 2, elem)>::type      \
  BOOST_PP_CAT (get_, BOOST_PP_TUPLE_ELEM (3, 0, elem)) (                     \
    const tools::named_pbase_value_type_helper <BOOST_PP_TUPLE_ELEM (3, 2, elem)>::type &def \
    ) const                                                                   \
  {                                                                           \
    idx_type index = BOOST_PP_TUPLE_ELEM (3, 0, elem);                        \
    return tools::named_pbase_value_type_helper <BOOST_PP_TUPLE_ELEM (3, 2, elem)>::get_d (this, index, def); \
  }                                                                           \
  bool                                                                        \
  BOOST_PP_CAT (check_, BOOST_PP_TUPLE_ELEM (3, 0, elem)) () const            \
  {                                                                           \
    idx_type index = BOOST_PP_TUPLE_ELEM (3, 0, elem);                        \
    return this->check_value (index);                                         \
  }
#define PARAMS_DECLARE_CLASS_EREG_I(r, data, i, elem)                         \
  ereg <data> (BOOST_PP_TUPLE_ELEM (3, 0, elem),                              \
    BOOST_PP_STRINGIZE (BOOST_PP_TUPLE_ELEM (3, 0, elem)),                    \
    BOOST_PP_TUPLE_ELEM (3, 1, elem),                                         \
    BOOST_PP_TUPLE_ELEM (3, 2, elem));

#define MAIN_PARAMS_DECLARE_CLASS(class_name, seq)                            \
  struct BS_API_PLUGIN class_name : public named_pbase                        \
  {                                                                           \
    PROP_BASE_IDX_DECL_BEGIN (class_name, named_pbase)                        \
    BOOST_PP_SEQ_FOR_EACH_I (PARAMS_DECLARE_CLASS_I, _, seq)                  \
    CLASS_NAME_TOTAL,                                                         \
    PROP_BASE_IDX_DECL_END                                                    \
    PBASE_ACCESS_MS (class_name);                                             \
    BOOST_PP_SEQ_FOR_EACH_I (PARAMS_DECLARE_CLASS_GETTERS_I, _, seq)          \
    ~class_name () {}                                                         \
    class_name (bs_type_ctor_param )                                          \
    : bs_refcounter (), named_pbase ()                                        \
    {                                                                         \
      this->resize (CLASS_NAME_TOTAL);                                        \
      BOOST_PP_SEQ_FOR_EACH_I (PARAMS_DECLARE_CLASS_EREG_I, class_name, seq)  \
    }                                                                         \
    class_name (const class_name &c)                                          \
    : bs_refcounter (c), named_pbase ()                                       \
    {                                                                         \
      this->resize (CLASS_NAME_TOTAL);                                        \
      BOOST_PP_SEQ_FOR_EACH_I (PARAMS_DECLARE_CLASS_EREG_I, class_name, seq)  \
    }                                                                         \
    PARAMS_DECL_BLUE_SKY (class_name, named_pbase, BOOST_PP_STRINGIZE (class_name))  \
  };

#define MAIN_PARAMS(seq)                                                      \
  MAIN_PARAMS_DECLARE_CLASS (main_params_class, seq);                         \
  smart_ptr <main_params_class, true>     main_params_;                       \
  smart_ptr <named_pbase, true>                                               \
  main_params () const                                                        \
  {                                                                           \
    return main_params_;                                                      \
  }

#define NEXT_LINE_PARAMS(seq)                                                 \
  MAIN_PARAMS_DECLARE_CLASS (next_line_params_class, seq);                    \
  typedef smart_ptr <next_line_params_class, true>  sp_next_line_params_t;    \
  shared_vector <sp_next_line_params_t>      next_line_params_;               \
  smart_ptr <named_pbase, true>                                               \
  add_next_line_params ()                                                     \
  {                                                                           \
    sp_next_line_params_t p = BS_KERNEL.create_object (next_line_params_class::bs_type ()); \
    if (!p)                                                                   \
      {                                                                       \
        bs_throw_exception ("Can't create next_line_params");                 \
      }                                                                       \
    next_line_params_.push_back (p);                                          \
    return p;                                                                 \
  }                                                                           \
  bool is_multi_line () const { return true; }



#endif // #ifndef BS_BOS_CORE_WELL_EVENT_PARAMS_DECL_H_
