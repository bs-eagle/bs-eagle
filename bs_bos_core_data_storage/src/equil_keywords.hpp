#ifndef EQUIL_KEYWORDS_HPP_ed5cc068_7247_11e0_93af_cfbf98c3a73f
#define EQUIL_KEYWORDS_HPP_ed5cc068_7247_11e0_93af_cfbf98c3a73f
/**
 *       \file  equil_keywords.hpp
 *      \brief  Keywords for EQUIL model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  29.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "keyword_info_base.h"

namespace blue_sky 
{
  class BS_API_PLUGIN equil_keywords : public keyword_info_base
  {
  public:
    BLUE_SKY_TYPE_DECL (equil_keywords);

    virtual void
    register_keywords (sp_objbase &keyword_manager, std::string provider) const;
  };
}

#endif //
