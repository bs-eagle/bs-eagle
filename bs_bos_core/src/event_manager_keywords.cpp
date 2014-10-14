/**
 *       \file  event_manager_keywords.cpp
 *      \brief  keywords like DATE, DATES, TSTEP, TSTEPS
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  05.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "bs_common.h"
#include BS_FORCE_PLUGIN_IMPORT()
#include "bos_reader_iface.h"
#include "event_manager_iface.hpp"
#include "event_manager.h"
#include BS_STOP_PLUGIN_IMPORT()
#include "event_manager_keywords.hpp"

namespace blue_sky 
{
  void
  DATE (std::string const &keyword, keyword_params &params)
  {
    BS_SP (bos_reader_iface) reader = params.hdm->get_reader ();
    BS_SP (event_manager_iface) em = params.hdm->get_event_manager ();

    char buf[4096] = {0};
    if (reader->read_line (buf, 4096) <= 0)
      {
        bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
          % reader->get_prefix ().c_str () % keyword);
      }

    event_manager::date_t date;
    if (reader->get_dt ()->cstr2d (buf, date))
      {
        bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
          % reader->get_prefix ().c_str () % keyword);
      }

    em->set_current_date (date);
    BOSOUT (section::read_data, level::medium) << keyword << " " << date << bs_end;
  }

  void
  TSTEP (std::string const &keyword, keyword_params &params)
  {
    BS_SP (bos_reader_iface) reader = params.hdm->get_reader ();
    BS_SP (event_manager_iface) em = params.hdm->get_event_manager ();

    boost::array <t_float, 4096> tmp;
    t_long len = reader->read_fp_array (keyword.c_str (), &tmp[0], 4096);
    if (len < 1)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix () % keyword);
      }

    // dtmp[i] <- days
    // dtmp[i]*24*60*60 <-seconds
    // dtmp[i]*24*60*60*1000 <-millisecs
    // FIXME: may be we can ommit unused dates and do set_current_date once after loop
    for (t_long i = 0; i < len; ++i)
      {
        event_manager::date_t date = em->get_current_date ();
        //date += boost::posix_time::millisec (tmp[i] * 24 * 60 * 60 * 1000);
        date += tmp[i];
        em->set_current_date (date);
      }

    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  event_manager_keywords::event_manager_keywords (bs_type_ctor_param)
  {
  }

  event_manager_keywords::event_manager_keywords (event_manager_keywords const &src)
  : bs_refcounter (src), keyword_info_base (src)
  {
  }

  void
  event_manager_keywords::register_keywords (sp_objbase &km, std::string /*provider*/) const
  {
    BS_SP (keyword_manager_iface) keyword_manager (km, bs_dynamic_cast ());
    BS_ASSERT (keyword_manager);

    keyword_manager->register_keyword ("DATE",    keyword_handler (DATE, 0));
    keyword_manager->register_keyword ("DATES",   keyword_handler (DATE, 0));
    keyword_manager->register_keyword ("TSTEP",   keyword_handler (TSTEP, 0));
    keyword_manager->register_keyword ("TSTEPS",  keyword_handler (TSTEP, 0));

  }

  BLUE_SKY_TYPE_STD_CREATE (event_manager_keywords);
  BLUE_SKY_TYPE_STD_COPY (event_manager_keywords);
  BLUE_SKY_TYPE_IMPL (event_manager_keywords, keyword_info_base, "Event Manager Keywords", "event_manager_keywords", "event_manager_keywords");
}
