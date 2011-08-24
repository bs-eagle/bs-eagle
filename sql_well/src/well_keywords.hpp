#ifndef WELL_KEYWORDS_H_
#define WELL_KEYWORDS_H_
/**
 *       \file  well_keywords.hpp
 *      \brief  keywords for wells
 *     \author  Mark Khait
 *       \date  05.08.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "keyword_info_base.h"

namespace blue_sky
{
  class BS_API_PLUGIN well_keywords : public keyword_info_base
  {
  public:
    typedef well_keywords             this_t;
    
  public:
    BLUE_SKY_TYPE_DECL (well_keywords);
        

    virtual void
    register_keywords (sp_objbase &keyword_manager, std::string provider) const;
    
    static void CSV_SHEDULE_reactor (std::string const &keyword, keyword_params &params);
  };
}

#endif //WELL_KEYWORDS_H_
