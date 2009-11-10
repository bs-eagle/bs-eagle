/*!
  \file err_num_def.h
  \brief Base error codes definitions
*/

#ifndef __ERROR_DEF
#define __ERROR_DEF

#define YS_SUCCESS                              (0)
#define YS_NO_MEMORY                            (-1)
#define YS_NO_ARGUMENTS                         (-2)
#define YS_BAD_INPUT_POINTER                    (-3)
#define YS_STACK_IS_EMPTY                       (-4)
#define YS_STACK_IS_FULL                        (-5)
#define YS_ERROR_POINTER_IS_ZERO                (-6)
#define YS_BAD_VALUE                            (-7)
#define YS_VALUES_NOT_EQUALS                    (-8)
#define YS_PVT_TABLES_IS_ZERO                   (-9)
#define YS_BOSATS_RETURN_ERROR                  (-10)
#define YS_CANNOT_ADD_WELL                      (-11)
#define YS_MINIMAL_PAV                          (-12)
#define YS_UNKNOWN_WELL_STATUS                  (-13)
#define YS_ZCOMP_RETURN_ERROR                   (-14)
#define YS_CANNOT_ADD_CONNECTION                (-15)
#define YS_NO_OPEN_FILE                         (-16)
#define YS_NOTHING_TO_WRITE                     (-17)
#define YS_CANNOT_OPEN_FILE                     (-18)
#define YS_SOIL_IS_OFR                          (-19)
#define YS_SWAT_IS_OFR                          (-20)
#define YS_SWOF_IS_OFR                          (-21)
#define YS_SGOF_IS_OFR                          (-22)
#define YS_VALUE_HAVE_NOT_INIT                  (-23)
#define YS_NOTHING_TO_READ                      (-25)
#define YS_CANNOT_READ_FROM_FILE                (-26)
#define YS_ORDER_ERROR                          (-27)
#define YS_ERROR_IN_WELL_NAME                   (-28)
#define YS_INDEX_OUT_OF_RANGE                   (-29)
#define YS_CALC_ERR                             (-30)
#define YS_CANNOT_ADD_GROUP                     (-31)
#define YS_CANNOT_FIND_CONNECTION               (-32)
#define YS_STREAMLINE_CALC                      (-33)
#define YS_EOF                                  (-34)

#define YS_SWOF_IS_NONE_INCREASING              (-35)
#define YS_SGOF_IS_NONE_INCREASING              (-36)
#define YS_SWFN_IS_NONE_INCREASING              (-37)
#define YS_SGFN_IS_NONE_INCREASING              (-38)
#define YS_SOF2_IS_NONE_INCREASING              (-39)
#define YS_SOF3_IS_NONE_INCREASING              (-40)
#define YS_SOF32D_IS_NONE_INCREASING            (-41)

// for UfaSolverisub
#define YS_NO_TIME_STEPS                        (-101)
#define YS_UNKNOWN_UNITS                        (-102)
#define YS_CURRENTLY_NOT_SUPPORTED_UNITS        (-103)

#define YS_SWFN_IS_OFR                          (-104)
#define YS_SGFN_IS_OFR                          (-105)
#define YS_SOF3_IS_OFR                          (-106)
#define YS_SOF2_IS_OFR                          (-107)
#define YS_SOF32D_IS_OFR                        (-108)

//! small chops on time step
#define YS_SMALL_TIME_STEP_CHOP                 (256)

// main definition
#define COMP_EPSILON 1e-15
#define CALC_EPSILON 1e-8

#ifdef _SGM_KERNEL_
#define UFASOLVER_VERSION_STRING ("SGM Kernel version 0.9.9")
#define UFASOLVER_SIGNATURE (/* GET_TEXT */"SGM Kernel")
#else // !_SGM_KERNEL_
#ifdef _SGM_KERNEL_DEMO_
#define UFASOLVER_VERSION_STRING ("SGM Kernel version 0.9.9")
#define UFASOLVER_SIGNATURE (/* GET_TEXT */"SGM Kernel (DEMO VERSION)")
#else // !_SGM_KERNEL_DEMO_
#ifndef _UFASOLVER_DEMO_VERSION_
#define UFASOLVER_VERSION_STRING ("UfaSolver version 0.9.8")
#define UFASOLVER_SIGNATURE (/*GET_TEXT*/"UfaSolver")
#else // _UFASOLVER_DEMO_VERSION_
#define UFASOLVER_VERSION_STRING ("UfaSolver version 0.9.8")
#define UFASOLVER_SIGNATURE (/*GET_TEXT*/"UfaSolver (DEMO VERSION)")
#endif // _UFASOLVER_DEMO_VERSION_
#endif // _SGM_KERNEL_DEMO_
#endif // _SGM_KERNEL_

#endif // __ERROR_DEF
