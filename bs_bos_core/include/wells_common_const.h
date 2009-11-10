/**
 * \file wells_common_const.h
 * \brief wells common consts holder
 * \author Sergey Miryanov
 * \date 07.07.2008
 * */
#ifndef BS_WELLS_COMMON_CONST_H_
#define BS_WELLS_COMMON_CONST_H_

namespace blue_sky
  {

#ifndef MIN_ZERO_DIFF
#define MIN_ZERO_DIFF                   1.0e-14
#endif

#ifndef DIFF_EPSILON
#define DIFF_EPSILON                    1.0e-12
#endif

#ifndef PI
#define PI                              3.141592653589793238463
#endif

#ifndef MIN_FRACTURE_MULT_CORR
#define MIN_FRACTURE_MULT_CORR          0.001
#endif

}


#endif  // #ifndef BS_WELLS_COMMON_CONST_H_
