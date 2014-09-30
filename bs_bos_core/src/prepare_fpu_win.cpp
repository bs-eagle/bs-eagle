/**
 *       \file  prepare_fpu_win.cpp
 *      \brief  Prepares FPU to work (enables exceptions, etc), windows version
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  10.09.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

/**
 * \file prepare_fpu_win.cpp
 * \brief prepare fpu to work (enable exception, etc)
 * \author Sergey Miryanov
 * \date 10.09.2009
 * */
#include "prepare_fpu.h"

#ifndef UNIX
#include <float.h>

namespace blue_sky {
namespace tools {

  void
  prepare_fpu::enable_exceptions ()
  {
    _clearfp ();

    unsigned int cw = _controlfp (0, 0);
    cw &= ~(EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID);
    unsigned int original = _controlfp (cw, MCW_EM);

    // to restore call
    // _controlfp (original, MCW_EM)
  }

} // namespace tools
} // namespace blue_sky


#endif // #ifndef UNIX

