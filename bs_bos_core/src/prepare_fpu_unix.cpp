/**
 *       \file  prepare_fpu_unix.cpp
 *      \brief  Prepares FPU to work (enables exceptions, etc), linux version
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  10.09.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "prepare_fpu.h"

#ifdef UNIX
#include <fenv.h>
#include <signal.h>

#include "bs_kernel_tools.h"

namespace blue_sky {
namespace tools {

  namespace detail {

    void 
    set_fpu (unsigned int mode)
    {
      asm ("fldcw %0" : : "m" (*&mode));
    }

    void
    fpe_handler (int sig_num)
    {
      bs_throw_exception ("FPU exception occured");
      //  << kernel_tools::get_backtrace (72)
      //  << bs_end;

      //signal (SIGFPE, SIG_DFL);
      //raise (sig_num);
    }

  } // namespace detail

  void
  prepare_fpu::enable_exceptions ()
  {
    //detail::set_fpu (0x27F);
    //feenableexcept (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#ifdef BS_BOS_CORE_USE_FPE_HANDLER
    //signal (SIGFPE, detail::fpe_handler);
#endif
  }

} // namespace tools
} // namespace blue_sky

#endif // #ifndef UNIX

