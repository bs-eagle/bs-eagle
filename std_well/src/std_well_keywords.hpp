#ifndef STD_WELL_KEYWORDS_HPP_2013e144_76fd_11e0_9429_77ff7c099e25
#define STD_WELL_KEYWORDS_HPP_2013e144_76fd_11e0_9429_77ff7c099e25
/**
 *       \file  std_well_keywords.hpp
 *      \brief  keywords for STD WELL
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  05.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "keyword_info_base.h"

namespace blue_sky 
{
  class BS_API_PLUGIN std_well_keywords : public keyword_info_base
  {
  public:
    BLUE_SKY_TYPE_DECL (std_well_keywords);

    virtual void
    register_keywords (sp_objbase &keyword_manager, std::string provider) const;
  };
}

#endif 

