/**
 *       \file  scal_keywords.cpp
 *      \brief  keywords for SCAL
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  03.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "bs_scal_stdafx.h"
#include "scal_keywords.h"
#include "keyword_manager_iface.h"
#include "data_class.h"
#include "read_class.h"

namespace blue_sky 
{

void 
SWOF (const std::string &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->swof.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SWOF table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
        size_t len = 0;
        if ((len = reader->read_table (keyword, (*idata->swof[i].main_data_), 4)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) << "scal_region=" << n_scal_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
void 
SGOF (const std::string &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sgof.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SGOF table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
        size_t len = 0;
        if ((len = reader->read_table (keyword, (*idata->sgof[i].main_data_), 4)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) << "scal_region=" << n_scal_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  
void 
SWFN (const std::string &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->swfn.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SWFN table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
        size_t len = 0;
        if ((len = reader->read_table (keyword, (*idata->swfn[i].main_data_), 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) << "scal_region=" << n_scal_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


void 
SGFN (const std::string &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sgfn.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SGFN table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
        size_t len = 0;
        if ((len = reader->read_table (keyword, (*idata->sgfn[i].main_data_), 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) << "scal_region=" << n_scal_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  
void 
SOF3 (const std::string &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sof3.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SOF3 table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
        size_t len = 0;
        if ((len = reader->read_table (keyword, (*idata->sof3[i].main_data_), 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) << "scal_region=" << n_scal_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }
  
  
void 
SOF2 (const std::string &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sof2.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SOF3 table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
        size_t len = 0;
        if ((len = reader->read_table (keyword, (*idata->sof2[i].main_data_), 2)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) << "scal_region=" << n_scal_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  scal_keywords::scal_keywords (bs_type_ctor_param)
  {
  }

  scal_keywords::scal_keywords (const scal_keywords &src)
  : bs_refcounter (src), keyword_info_base (src)
  {
  }

  void
  scal_keywords::register_keywords (sp_objbase &km, std::string provider) const
  {
    BS_SP (keyword_manager_iface) keyword_manager (km, bs_dynamic_cast ());
    BS_ASSERT (keyword_manager);

    keyword_manager->register_keyword ("SWOF", keyword_handler (SWOF, 0));
    keyword_manager->register_keyword ("SGOF", keyword_handler (SGOF, 0));
    keyword_manager->register_keyword ("SWFN", keyword_handler (SWFN, 0));
    keyword_manager->register_keyword ("SGFN", keyword_handler (SGFN, 0));
    keyword_manager->register_keyword ("SOF2", keyword_handler (SOF2, 0));
    keyword_manager->register_keyword ("SOF3", keyword_handler (SOF3, 0));
  }


  BLUE_SKY_TYPE_STD_CREATE (scal_keywords);
  BLUE_SKY_TYPE_STD_COPY (scal_keywords);
  BLUE_SKY_TYPE_IMPL (scal_keywords, keyword_info_base, "SCAL Keywords", "scal_keywords", "scal_keywords");
}
