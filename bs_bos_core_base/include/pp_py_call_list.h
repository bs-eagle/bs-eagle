/**
 * \file pp_py_call_list.h
 * \brief
 * \author Sergey Miryanov
 * \date 26.06.2009
 * */
#ifndef BS_BOS_CORE_BASE_PP_PY_CALL_LIST_H_
#define BS_BOS_CORE_BASE_PP_PY_CALL_LIST_H_

#include "pp_param_list.h"

/**
 * */
#define PY_CALL_LIST(tuple_, size_) BOOST_PP_CAT (PY_CALL_LIST_, size_) (tuple_, size_, size_)
#define PY_CALL_LIST_0
#define PY_CALL_LIST_1(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_))
#define PY_CALL_LIST_2(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_1  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_3(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_2  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_4(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_3  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_5(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_4  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_6(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_5  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_7(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_6  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_8(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_7  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_9(tuple_, size_, rsize_)  pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_8  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_10(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_9  (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_11(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_10 (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_12(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_11 (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_13(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_12 (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_14(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_13 (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_15(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_14 (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_16(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_15 (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_17(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_16 (tuple_, size_, BOOST_PP_DEC (rsize_))
#define PY_CALL_LIST_18(tuple_, size_, rsize_) pass_arg_to_python (BS_PP_ELEM_NAME (A, size_, rsize_)), PY_CALL_LIST_17 (tuple_, size_, BOOST_PP_DEC (rsize_))

#endif // #ifndef BS_BOS_CORE_BASE_PP_PY_CALL_LIST_H_

