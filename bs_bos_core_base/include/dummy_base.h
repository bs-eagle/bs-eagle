/**
 * \file dummy_base.h
 * \brief dummy class for any base python wrappers
 * \author Miryanov Sergey
 * \date 08.05.2008
 */
#ifndef BS_DUMMY_BASE_H_
#define BS_DUMMY_BASE_H_

#include "bs_object_base.h"
#include "bs_link.h"

namespace blue_sky
  {

  /**
   * \brief dummy linear_solver base
   */
  class BS_API_PLUGIN dummy_base : public objbase
  {
  public:

    BLUE_SKY_TYPE_DECL (dummy_base);
  };


} // namespace blue_sky

#endif // #ifndef BS_DUMMY_BASE_H_
