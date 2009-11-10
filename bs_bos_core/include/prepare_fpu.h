/**
 * \file prepare_fpu.h
 * \brief prepare fpu to work (enable exception, etc)
 * \author Sergey Miryanov
 * \date 10.09.2009
 * */

namespace blue_sky {
namespace tools {

  struct prepare_fpu
  {
    static void
    enable_exceptions ();
  };


} // namespace tools
} // namespace blue_sky

