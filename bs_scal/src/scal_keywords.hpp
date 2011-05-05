#ifndef SCAL_KEYWORDS_H_97d0602c_70bd_11e0_9b92_43fd738d2557
#define SCAL_KEYWORDS_H_97d0602c_70bd_11e0_9b92_43fd738d2557
/**
 *       \file  scal_keywords.hpp
 *      \brief  keywords for SCAL
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  27.04.2011
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
  };
}

#endif //
