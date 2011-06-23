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

  keyword_info_base::keyword_info_base(bs_type_ctor_param)
    {

    }

  keyword_info_base::keyword_info_base(const keyword_info_base& src)
  : bs_refcounter (src), objbase (src)
    {
      // TODO: BUG:
      bs_throw_exception ("NOT IMPL YET");
      //*this = src;
    }

  BLUE_SKY_TYPE_IMPL (keyword_info_base, objbase, "keyword_info_base", "keyword_info_base", "keyword_info_base");

  BLUE_SKY_TYPE_STD_CREATE (keyword_info_base);
  BLUE_SKY_TYPE_STD_COPY (keyword_info_base);
}//ns bs

