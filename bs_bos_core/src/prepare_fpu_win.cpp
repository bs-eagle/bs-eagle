/**
 * \file prepare_fpu_win.cpp
 * \brief prepare fpu to work (enable exception, etc)
 * \author Sergey Miryanov
 * \date 10.09.2009
 * */
#include "stdafx.h"
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

