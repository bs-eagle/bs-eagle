/**
 *       \file  std_well_keywords.cpp
 *      \brief  keywords for STD WELL
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  05.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "std_well_stdafx.hpp"
#include "std_well_keywords.hpp"
#include "keyword_manager_iface.h"
#include "read_class.h"
#include "event_manager_iface.hpp"

namespace blue_sky
{
  void
  event_parser (std::string const &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (event_manager_iface) em = params.hdm->get_event_manager ();

    for (;;)
      {
        char buf[4096] = {0};
        t_long len = reader->read_line (buf, 4096);
        if (len == 0 || buf[0] == '/') 
          {
            em->end_event ();
            break;
          }

        em->process_event (em->get_current_date (), keyword, buf);
      }

  }

  std_well_keywords::std_well_keywords (bs_type_ctor_param)
  {
  }

  std_well_keywords::std_well_keywords (std_well_keywords const &src)
  : bs_refcounter (src), keyword_info_base (src)
  {
  }

  void
  std_well_keywords::register_keywords (sp_objbase &km, std::string provider) const
  {
    BS_SP (keyword_manager_iface) keyword_manager (km, bs_dynamic_cast ());
    BS_ASSERT (keyword_manager);

    keyword_manager->register_keyword ("WELSPECS",      keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WELLCON",       keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("COMPDAT",       keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WCONPROD",      keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WCONHIST",      keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WCONINJE",      keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WECON",         keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WECONINJ",      keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WEFAC",         keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WELTARG",       keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("WPIMULT",       keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("COMPENSATION",  keyword_handler (event_parser, 0));
    keyword_manager->register_keyword ("PERMFRAC",      keyword_handler (event_parser, 0));
  }

  BLUE_SKY_TYPE_STD_CREATE (std_well_keywords);
  BLUE_SKY_TYPE_STD_COPY (std_well_keywords);
  BLUE_SKY_TYPE_IMPL (std_well_keywords, keyword_info_base, "STD WELL Keywords", "STD_WELL", "STD_WELL");
}

