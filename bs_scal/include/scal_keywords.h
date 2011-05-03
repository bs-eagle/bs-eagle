#ifndef SCAL_KEYWORDS_H_599eaa62_70a3_11e0_b06f_475542e35bfb
#define SCAL_KEYWORDS_H_599eaa62_70a3_11e0_b06f_475542e35bfb
/**
 *       \file  scal_keywords.h
 *      \brief  keyword's handler for SCAL
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  03.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "keyword_info_base.h"

namespace blue_sky
{
  class BS_API_PLUGIN scal_keywords : public keyword_info_base
  {
  public:
    BLUE_SKY_TYPE_DECL (scal_keywords);

    virtual void
    register_keywords (sp_objbase &keyword_manager, std::string provider) const;

    //void
    //activate_keywords (BS_SP (keyword_manager_iface) keyword_manager);
  };
}

#endif //