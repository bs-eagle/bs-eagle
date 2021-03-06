/** 
 * \file python_method_wrapper.h
 * \brief
 * \author Sergey Miryanov
 * \date 22.05.2009
 * */
#ifndef BS_BOS_CORE_BASE_PYTHON_WRAPPER_METHOD_H_
#define BS_BOS_CORE_BASE_PYTHON_WRAPPER_METHOD_H_

#include "pass_arg_to_python.h"
#include "pp_py_call_list.h"
#include "throw_exception.h"


/** 
 * @brief wrap functions with out return type
 * 
 * @param name_         -- function name
 * @param ret_          -- should be void
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAPPER_METHOD(name_, ret_, len_, seq_)                                                       \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_))                                                                     \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        f (PY_CALL_LIST (seq_, len_));                                                                \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        wrapped_t::name_ (CALL_LIST (seq_, len_));                                                    \
      }                                                                                               \
  }                                                                                                   \
  ret_                                                                                                \
  BOOST_PP_CAT (wrapped_, name_) (PARAM_LIST (seq_, len_))                                            \
  {                                                                                                   \
    wrapped_t::name_ (CALL_LIST (seq_, len_));                                                        \
  }

/** 
 * @brief wrap constant functions with out return type
 * 
 * @param name_         -- function name
 * @param ret_          -- should be void
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAPPER_METHOD_CONST(name_, ret_, len_, seq_)                                                 \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_)) const                                                               \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        f (PY_CALL_LIST (seq_, len_));                                                                \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        wrapped_t::name_ (CALL_LIST (seq_, len_));                                                    \
      }                                                                                               \
  }                                                                                                   \
  ret_                                                                                                \
  BOOST_PP_CAT (wrapped_, name_) (PARAM_LIST (seq_, len_)) const                                      \
  {                                                                                                   \
    wrapped_t::name_ (CALL_LIST (seq_, len_));                                                        \
  }

/** 
 * @brief wrap functions
 * 
 * @param name_         -- function name
 * @param ret_          -- return type
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAPPER_METHOD_R(name_, ret_, len_, seq_)                                                     \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_))                                                                     \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        return f (PY_CALL_LIST (seq_, len_));                                                         \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        return wrapped_t::name_ (CALL_LIST (seq_, len_));                                             \
      }                                                                                               \
  }                                                                                                   \
  ret_                                                                                                \
  BOOST_PP_CAT (wrapped_, name_) (PARAM_LIST (seq_, len_))                                            \
  {                                                                                                   \
    return wrapped_t::name_ (CALL_LIST (seq_, len_));                                                 \
  }

/** 
 * @brief wrap constant functions
 * 
 * @param name_         -- function name
 * @param ret_          -- return type
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAPPER_METHOD_R_CONST(name_, ret_, len_, seq_)                                               \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_)) const                                                               \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        return f (PY_CALL_LIST (seq_, len_));                                                         \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        return wrapped_t::name_ (CALL_LIST (seq_, len_));                                             \
      }                                                                                               \
  }                                                                                                   \
  ret_                                                                                                \
  BOOST_PP_CAT (wrapped_, name_) (PARAM_LIST (seq_, len_)) const                                      \
  {                                                                                                   \
    return wrapped_t::name_ (CALL_LIST (seq_, len_));                                                 \
  }

/** 
 * @brief wrap pure virtual method with out return type
 * 
 * @param name_         -- function name
 * @param ret_          -- should be void
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAP_PURE_METHOD(name_, ret_, len_, seq_)                                                     \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_))                                                                     \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        f (PY_CALL_LIST (seq_, len_));                                                                \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        bs_throw_exception ("PURE_CALL");                                                             \
      }                                                                                               \
  }

/** 
 * @brief wrap pure virtual const method with out return type
 * 
 * @param name_         -- function name
 * @param ret_          -- should be void
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAP_PURE_METHOD_CONST(name_, ret_, len_, seq_)                                               \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_)) const                                                               \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        f (PY_CALL_LIST (seq_, len_));                                                                \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        bs_throw_exception ("PURE_CALL");                                                             \
      }                                                                                               \
  }

/** 
 * @brief wrap pure virtual method
 * 
 * @param name_         -- function name
 * @param ret_          -- return type
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAP_PURE_METHOD_R(name_, ret_, len_, seq_)                                                   \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_))                                                                     \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        return f (PY_CALL_LIST (seq_, len_));                                                         \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        bs_throw_exception ("PURE_CALL");                                                             \
      }                                                                                               \
  }

/** 
 * @brief wrap pure virtual const method
 * 
 * @param name_         -- function name
 * @param ret_          -- return type
 * @param len_          -- number of parameters
 * @param seq_          -- list of parameters
 * 
 * @return nothing
 */
#define WRAP_PURE_METHOD_R_CONST(name_, ret_, len_, seq_)                                             \
  ret_                                                                                                \
  name_ (PARAM_LIST (seq_, len_)) const                                                               \
  {                                                                                                   \
    using namespace blue_sky::python::tools;                                                          \
    if (boost::python::override f = this->get_override (#name_))                                      \
      {                                                                                               \
        return f (PY_CALL_LIST (seq_, len_));                                                         \
      }                                                                                               \
    else                                                                                              \
      {                                                                                               \
        bs_throw_exception ("PURE_CALL");                                                             \
      }                                                                                               \
  }
#endif // #ifndef BS_BOS_CORE_BASE_PYTHON_WRAPPER_METHOD_H_
