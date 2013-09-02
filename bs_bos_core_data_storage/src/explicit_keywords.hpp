#ifndef EXPLICIT_MODEL_HPP_60db5192_7254_11e0_80a3_7fc2a5b0f137
#define EXPLICIT_MODEL_HPP_60db5192_7254_11e0_80a3_7fc2a5b0f137
/**
 *       \file  explicit_keywords.hpp
 *      \brief  Keywords for EXPLICIT model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  29.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "keyword_info_base.h"

namespace blue_sky 
{
  class BS_API_PLUGIN explicit_keywords : public keyword_info_base
  {
  public:
    BLUE_SKY_TYPE_DECL (explicit_keywords);

    virtual void
    register_keywords (sp_objbase &keyword_manager, std::string provider) const;
  };
}

#endif //
