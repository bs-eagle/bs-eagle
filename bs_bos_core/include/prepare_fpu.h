/**
 *       \file  prepare_fpu.h
 *      \brief  Prepares FPU to work (enables exceptions, etc)
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  10.09.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
namespace blue_sky {
namespace tools {

  /**
   * \class prepare_fpu
   * \brief Prepares FPU to work
   * */
  struct prepare_fpu
  {
    /**
     * \brief  Enables exceptions
     * */
    static void
    enable_exceptions ();
  };


} // namespace tools
} // namespace blue_sky

