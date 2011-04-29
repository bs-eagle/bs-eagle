/**
 *       \file  equil_keywords.cpp
 *      \brief  Keywords for EQUIL model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  29.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "equil_keywords.hpp"
#include "keyword_manager_iface.h"

namespace blue_sky 
{
  equil_keywords::equil_keywords (bs_type_ctor_param)
  {
  }

  equil_keywords::equil_keywords (equil_keywords const &src)
  : bs_refcounter (src), keyword_info_base (src)
  {
  }

  namespace 
  {
    void
    EQUIL (std::string const &, keyword_params &params)
    {
    }
    void
    RSVD (std::string const &, keyword_params &params)
    {
    }
    void
    PBVD (std::string const &, keyword_params &params)
    {
    }

    void
    activate (std::string const &, keyword_params &params)
    {
      BS_SP (keyword_manager_iface) keyword_manager = params.hdm->get_keyword_manager ();
      BS_ASSERT (keyword_manager);

      keyword_manager->register_keyword ("EQUIL", keyword_manager_iface::keyword_handler (EQUIL, 0));
      keyword_manager->register_keyword ("RSVD", keyword_manager_iface::keyword_handler (RSVD, 0));
      keyword_manager->register_keyword ("PBVD", keyword_manager_iface::keyword_handler (PBVD, 0));
    }
  }

  void
  equil_keywords::register_keywords (sp_objbase &km, std::string provider) const
  {
    BS_SP (keyword_manager_iface) keyword_manager (km, bs_dynamic_cast ());
    BS_ASSERT (keyword_manager);

    if (provider == "")
      {
        provider = "EQUIL_MODEL";
        keyword_manager->register_keyword ("EQUIL_MODEL", keyword_handler (0, activate));
      }

    keyword_manager->register_supported_keyword ("EQUIL", provider);
    keyword_manager->register_supported_keyword ("RSVD", provider);
    keyword_manager->register_supported_keyword ("PBVD", provider);
  }

  BLUE_SKY_TYPE_STD_CREATE (equil_keywords);
  BLUE_SKY_TYPE_STD_COPY (equil_keywords);
  BLUE_SKY_TYPE_IMPL (equil_keywords, keyword_info_base, "Keywords for EQUIL model", "EQUIL_MODEL", "EQUIL_MODEL");
}
