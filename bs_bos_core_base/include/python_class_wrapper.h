/**
 * \file python_class_wrapper.h
 * \brief
 * \author Sergey Miryanov
 * \date 26.06.2009
 * */
#ifndef BS_BOS_CORE_BASE_PYTHON_CLASS_WRAPPER_H_
#define BS_BOS_CORE_BASE_PYTHON_CLASS_WRAPPER_H_

#include "make_me_happy.h"
#include "pp_param_list.h"

#define CLASS_WRAPPER_BASE_NAME(S)     BOOST_PP_CAT (S, _base)

#define CLASS_WRAPPER(class_name_, wrapper_name_)                                                                                   \
    struct BS_API_PLUGIN CLASS_WRAPPER_BASE_NAME (wrapper_name_): public class_name_                                                \
    {                                                                                                                               \
      MAKE_ME_HAPPY (CLASS_WRAPPER_BASE_NAME (wrapper_name_), class_name_, BOOST_PP_STRINGIZE (CLASS_WRAPPER_BASE_NAME (wrapper_name_)));  \
    };                                                                                                                              \
    class BS_API_PLUGIN wrapper_name_ : public CLASS_WRAPPER_BASE_NAME (wrapper_name_), public boost::python::wrapper <CLASS_WRAPPER_BASE_NAME (wrapper_name_)>

#define CLASS_WRAPPER_T(type_count_, type_list_, class_name_, wrapper_name_)                                                        \
  template <TYPE_LIST_DECL (type_list_, type_count_)>                                                                               \
  struct BS_API_PLUGIN CLASS_WRAPPER_BASE_NAME (wrapper_name_): public class_name_ <TYPE_LIST (type_list_, type_count_)>            \
  {                                                                                                                                 \
    MAKE_ME_HAPPY (CLASS_WRAPPER_BASE_NAME (wrapper_name_), class_name_ <TYPE_LIST (type_list_, type_count_)>, BOOST_PP_STRINGIZE (CLASS_WRAPPER_BASE_NAME (wrapper_name_)));  \
  };                                                                                                                                \
  template <TYPE_LIST_DECL (type_list_, type_count_)>                                                                               \
  class BS_API_PLUGIN wrapper_name_ : public CLASS_WRAPPER_BASE_NAME (wrapper_name_) <TYPE_LIST (type_list_, type_count_)>, public boost::python::wrapper <CLASS_WRAPPER_BASE_NAME (wrapper_name_) <TYPE_LIST (type_list_, type_count_)> >

#define CLASS_WRAPPER_DECL_T(type_count_, type_list_, wrapper_name_)                                                                \
  MAKE_ME_HAPPY (wrapper_name_, CLASS_WRAPPER_BASE_NAME (wrapper_name_) <TYPE_LIST (type_list_, type_count_)>, BOOST_PP_STRINGIZE (wrapper_name_));

#define CLASS_WRAPPER_DECL(wrapper_name_)                                                                                           \
  MAKE_ME_HAPPY (wrapper_name_, CLASS_WRAPPER_BASE_NAME (wrapper_name_), BOOST_PP_STRINGIZE (wrapper_name_));

#define STRATEGY_CLASS_WRAPPER(class_name_, wrapper_name_)                                                                          \
  CLASS_WRAPPER_T (1, (strategy_t), class_name_, wrapper_name_)

#define STRATEGY_CLASS_WRAPPER_DECL(wrapper_name_)                                                                                  \
  CLASS_WRAPPER_DECL_T (1, (strategy_t), wrapper_name_)

#define MATRIX_IFACE_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (2, (fp_vector_type, i_vector_type), class_name_, wrapper_name_)

#define MATRIX_IFACE_WRAPPER_DECL(wrapper_name_)                                                                              \
  CLASS_WRAPPER_DECL_T (2, (fp_vector_type, i_vector_type), wrapper_name_)

#define BCSR_MATRIX_IFACE_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), class_name_, wrapper_name_)

#define BDIAG_MATRIX_IFACE_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), class_name_, wrapper_name_)

#define JAC_MATRIX_IFACE_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), class_name_, wrapper_name_)

#define BCSR_MATRIX_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), class_name_, wrapper_name_)

#define BDIAG_MATRIX_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), class_name_, wrapper_name_)

#define BCSR_MATRIX_IFACE_WRAPPER_DECL(wrapper_name_)                                                                              \
  CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), wrapper_name_)

#define BDIAG_MATRIX_IFACE_WRAPPER_DECL(wrapper_name_)                                                                              \
  CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), wrapper_name_)

#define JAC_MATRIX_IFACE_WRAPPER_DECL(wrapper_name_)                                                                              \
  CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), wrapper_name_)

#define BCSR_MATRIX_WRAPPER_DECL(wrapper_name_)                                                                              \
  CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), wrapper_name_)

#define MATRIX_CLASS_WRAPPER(class_name_, wrapper_name_)                                                                            \
  CLASS_WRAPPER_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), class_name_, wrapper_name_)

#define MATRIX_CLASS_WRAPPER_DECL(wrapper_name_)                                                                                    \
  CLASS_WRAPPER_DECL_T (3, (fp_vector_type, i_vector_type, fp_storage_vector_type), wrapper_name_)


#define PROP_IFACE_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (4, (fp_type_t, i_type_t, s_type_t, b_type_t), class_name_, wrapper_name_)

#define PROP_IFACE_WRAPPER_DECL(wrapper_name_)                                                                              \
  CLASS_WRAPPER_DECL_T (4, (fp_type_t, i_type_t, s_type_t, b_type_t), wrapper_name_)

#define LSOLVER_IFACE_WRAPPER(class_name_, wrapper_name_)                                                                      \
  CLASS_WRAPPER_T (1, (strategy_t), class_name_, wrapper_name_)

#define LSOLVER_IFACE_WRAPPER_DECL(wrapper_name_)                                                                              \
  CLASS_WRAPPER_DECL_T (1, (strategy_t), wrapper_name_)
#endif // #ifndef BS_BOS_CORE_BASE_PYTHON_CLASS_WRAPPER_H_

