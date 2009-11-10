/**
* @file keywords.cpp
* @brief keywords info base class
* @author Mark Khait
* @date 2009-08-03
* */
#include "bs_bos_core_base_stdafx.h"
#include "keyword_info_base.h"
#include "strategies.h"

namespace blue_sky {

  template<class strategy_t>
  keyword_info_base<strategy_t>::keyword_info_base(bs_type_ctor_param)
    {

    }

  template<class strategy_t>
  keyword_info_base<strategy_t>::keyword_info_base(const keyword_info_base<strategy_t>& src)
  : bs_refcounter (src), objbase (src)
    {
      // TODO: BUG:
      bs_throw_exception ("NOT IMPL YET");
      //*this = src;
    }

  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (keyword_info_base<base_strategy_fi>) , 1, (objbase), "keyword_info_base_fi", "Keyword info base (virtual) class", "Keyword info base (virtual) class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (keyword_info_base<base_strategy_di>) , 1, (objbase), "keyword_info_base_di", "Keyword info base (virtual)  class", "Keyword info base (virtual)  class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (keyword_info_base<base_strategy_mixi>) , 1, (objbase), "keyword_info_base_mixi", "Keyword info base (virtual)  class", "Keyword info base (virtual)  class", false);

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(keyword_info_base, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(keyword_info_base, (class));
}//ns bs

