#ifndef PVT_KEYWORDS_H_599eaa62_70a3_11e0_b06f_475542e35bfb
#define PVT_KEYWORDS_H_599eaa62_70a3_11e0_b06f_475542e35bfb
/**
 *       \file  pvt_keywords.hpp
 *      \brief  keyword for PVT
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  27.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "keyword_info_base.h"

namespace blue_sky
{
  class BS_API_PLUGIN pvt_keywords : public keyword_info_base
  {
  public:
    BLUE_SKY_TYPE_DECL (pvt_keywords);

    virtual void
    register_keywords (sp_objbase &keyword_manager, std::string provider) const;

    //void
    //activate_keywords (BS_SP (keyword_manager_iface) keyword_manager);
  };
}

#endif //
