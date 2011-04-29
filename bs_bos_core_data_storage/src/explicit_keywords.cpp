/**
 *       \file  explicit_keywords.cpp
 *      \brief  Keywords for EXPLICIT model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  29.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "explicit_keywords.hpp"
#include "keyword_manager_iface.h"

namespace blue_sky 
{
  explicit_keywords::explicit_keywords (bs_type_ctor_param)
  {
  }

  explicit_keywords::explicit_keywords (explicit_keywords const &src)
  : bs_refcounter (src), keyword_info_base (src)
  {
  }

  namespace 
  {
    void
    PRVD (std::string const &, keyword_params &params)
    {
    }

    void
    activate (std::string const &, keyword_params &params)
    {
      BS_SP (keyword_manager_iface) keyword_manager = params.hdm->get_keyword_manager ();
      BS_ASSERT (keyword_manager);

      // FIXME: npy_intp
      npy_intp dimens[] = {1, 0, 1, 0, 1, 0};
      keyword_manager->register_fp_pool_keyword ("PRESSURE", dimens, 0, 0);
      keyword_manager->register_fp_pool_keyword ("SWAT",    dimens, 0, 0);
      keyword_manager->register_fp_pool_keyword ("SGAS",    dimens, 0, 0);
      keyword_manager->register_fp_pool_keyword ("SOIL",    dimens, 0, 0);
      keyword_manager->register_fp_pool_keyword ("RS",      dimens, 0, 0);
      keyword_manager->register_fp_pool_keyword ("PBUB",    dimens, 0, 0);

      // FIXME: if PRVD present in model PRESSURE shouldn't be specified
      keyword_manager->register_keyword ("PRVD", keyword_manager_iface::keyword_handler (PRVD, 0));
    }
  }

  void
  explicit_keywords::register_keywords (sp_objbase &km, std::string provider) const
  {
    BS_SP (keyword_manager_iface) keyword_manager (km, bs_dynamic_cast ());
    BS_ASSERT (keyword_manager);

    if (provider == "")
      {
        provider = "EXPLICIT_MODEL";
        keyword_manager->register_keyword ("EXPLICIT_MODEL", keyword_handler (0, activate));
      }

    keyword_manager->register_supported_keyword ("PRESSURE", provider);
    keyword_manager->register_supported_keyword ("PRVD", provider);
    keyword_manager->register_supported_keyword ("SWAT", provider);
    keyword_manager->register_supported_keyword ("SGAS", provider);
    keyword_manager->register_supported_keyword ("SOIL", provider);
    keyword_manager->register_supported_keyword ("RS", provider);
    keyword_manager->register_supported_keyword ("PBUB", provider);
  }

  BLUE_SKY_TYPE_STD_CREATE (explicit_keywords);
  BLUE_SKY_TYPE_STD_COPY (explicit_keywords);
  BLUE_SKY_TYPE_IMPL (explicit_keywords, keyword_info_base, "Keywords for EXPLCIIT model", "EXPLICIT_MODEL", "EXPLICIT_MODEL");

}
