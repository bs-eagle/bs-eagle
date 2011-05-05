#ifndef EVENT_MANAGER_KEYWORDS_HPP_1a5c5566_7715_11e0_bc6b_13290d8e0e0b
#define EVENT_MANAGER_KEYWORDS_HPP_1a5c5566_7715_11e0_bc6b_13290d8e0e0b
/**
 *       \file  event_manager_keywords.hpp
 *      \brief  keywords like DATE, DATES, TSTEP, TSTEPS
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  05.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "keyword_info_base.h"

namespace blue_sky 
{
  class BS_API_PLUGIN event_manager_keywords : public keyword_info_base
  {
  public:
    BLUE_SKY_TYPE_DECL (event_manager_keywords);

    virtual void
    register_keywords (sp_objbase &keyword_manager, std::string provider) const;
  };
}


#endif

