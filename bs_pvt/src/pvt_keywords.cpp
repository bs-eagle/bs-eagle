/**
 *       \file  pvt_keywords.cpp
 *      \brief  keywords for PVT
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  27.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "bs_pvt_stdafx.h"
#include "pvt_keywords.hpp"
#include "keyword_manager_iface.h"
#include "data_class.h"
#include "read_class.h"
#include "pvt_3p_iface.h"
#include "table_iface.h"

namespace blue_sky 
{
  void
  DENSITY (std::string const &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (pvt_3p_iface) pvt = params.hdm->get_pvt ();

    t_long n_pvt_region = params.hdm->get_prop()->get_i ("pvt_region");
    spv_float sp_density = pvt->get_density();

    sp_density->resize (3 * n_pvt_region);
    for (t_long i = 0; i < n_pvt_region; ++i)
      {
        if (reader->read_array (keyword, *sp_density, 3 * i, 3) != 3)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
      }

    
    pvt->set_density_to_pvt_internal();
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  void
  PVTO (std::string const &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (pvt_3p_iface) pvt = params.hdm->get_pvt ();
    BS_SP(table_iface) tbl;
    
    stdv_float dbuf;
    t_long n_pvt_region = params.hdm->get_prop()->get_i("pvt_region");
    dbuf.resize(4);
    /*   
    if (!idata->pvto.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for oil has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
    */

    // Read table for each of region
    for (t_long i = 0; i < n_pvt_region; i++)
      {
        tbl = pvt->get_table (i, FI_PHASE_OIL);
        tbl->init(0, 4);
        //t_long lj = 0;
        char buf[CHAR_BUF_LEN] = {0};
        if (reader->read_line (buf, CHAR_BUF_LEN) <= 0)
          {
            bs_throw_exception (boost::format ("Error in %s: can't read argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        if (sscanf (buf, "%lf%lf%lf%lf", &dbuf[0],
                    &dbuf[1], &dbuf[2],
                    &dbuf[3]) != 4)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        tbl->push_back(dbuf);
        for (;;)
          {
            if (reader->read_line (buf, CHAR_BUF_LEN) <= 0)
              {
                bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
                  % reader->get_prefix () % keyword);
              }
            if (buf[0] == '/')
              {
                break;
              }

            size_t len = 0;
            if ((len = sscanf (buf, "%lf%lf%lf%lf", &dbuf[0],
                               &dbuf[1], &dbuf[2],
                               &dbuf[3])) != 4)
              {
                if (len == 3)
                  {
                    // 3 values table
                    dbuf[3] = dbuf[2];
                    dbuf[2] = dbuf[1];
                    dbuf[1] = dbuf[0];
                    dbuf[0] = -999.0; 
                  }
                else
                  {
                    bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
                      % reader->get_prefix () % keyword);
                  }
              }
            tbl->push_back(dbuf);
          }

        BOSOUT (section::read_data, level::medium) << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }
  
  
  void
  PVDO (std::string const &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (pvt_3p_iface) pvt = params.hdm->get_pvt ();
    BS_SP(table_iface) tbl;
    int n_cols = 3;
    
    stdv_float dbuf;
    char buf[CHAR_BUF_LEN] = {0};
    t_long n_pvt_region = params.hdm->get_prop()->get_i("pvt_region");
    
    dbuf.resize(n_cols);
    
    // Read table for each of region
    for (t_int i = 0; i < n_pvt_region; i++)
      {
        tbl = pvt->get_table (i, FI_PHASE_OIL);
        tbl->init(0, n_cols);
        
        
        for (;;)
          {
            if (reader->read_line (buf, CHAR_BUF_LEN) <= 0)
              {
                bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
                  % reader->get_prefix () % keyword);
              }
            if (buf[0] == '/')
              {
                break;
              }

            size_t len = 0;
            if ((len = sscanf (buf, "%lf%lf%lf", &dbuf[0],
                               &dbuf[1], &dbuf[2])) != n_cols)
                  {
                    bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
                      % reader->get_prefix () % keyword);
                  }
            tbl->push_back(dbuf);
          }
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }
  
  void
  PVTW (std::string const &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (pvt_3p_iface) pvt = params.hdm->get_pvt ();
    BS_SP(table_iface) tbl;
    int n_rows = 0, n_cols = 4;
    
    stdv_float dbuf;
    char buf[CHAR_BUF_LEN] = {0};
    t_long n_pvt_region = params.hdm->get_prop()->get_i("pvt_region");
    
    dbuf.resize(n_cols);
    
    // Read table for each of region
    for (t_int i = 0; i < n_pvt_region; i++)
      {
        tbl = pvt->get_table (i, FI_PHASE_WATER);
        tbl->init(0, n_cols);
        n_rows = 0;
        
        for (;;n_rows++)
          {
            if (reader->read_line (buf, CHAR_BUF_LEN) <= 0)
              {
                bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
                  % reader->get_prefix () % keyword);
              }
            if (buf[0] == '/')
              {
                break;
              }
            if (n_rows > 1)
              {
                bs_throw_exception (boost::format ("Error in %s: you can specify only 4 arguments for keyword %s")
                  % reader->get_prefix () % keyword);
              }

            if (sscanf (buf, "%lf%lf%lf%lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3]) != n_cols)
                  {
                    bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
                      % reader->get_prefix () % keyword);
                  }
            tbl->push_back(dbuf);
          }
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  void
  PVDG (std::string const &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (pvt_3p_iface) pvt = params.hdm->get_pvt ();
    BS_SP(table_iface) tbl;
    int n_cols = 3;
    
    stdv_float dbuf;
    char buf[CHAR_BUF_LEN] = {0};
    t_long n_pvt_region = params.hdm->get_prop()->get_i("pvt_region");
    
    dbuf.resize(n_cols);
    
    // Read table for each of region
    for (t_int i = 0; i < n_pvt_region; i++)
      {
        tbl = pvt->get_table (i, FI_PHASE_GAS);
        tbl->init(0, n_cols);
        
        
        for (;;)
          {
            if (reader->read_line (buf, CHAR_BUF_LEN) <= 0)
              {
                bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
                  % reader->get_prefix () % keyword);
              }
            if (buf[0] == '/')
              {
                break;
              }

            if (sscanf (buf, "%lf%lf%lf", &dbuf[0], &dbuf[1], &dbuf[2]) != n_cols)
                  {
                    bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
                      % reader->get_prefix () % keyword);
                  }
            tbl->push_back(dbuf);
          }
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }
  
  void
  ROCK (std::string const &keyword, keyword_params &params)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();

    t_float *p_ref = idata->p_ref->data ();
    t_float*rock = idata->rock->data ();
    t_long n_pvt_region = idata->props->get_i("pvt_region");
    
    for (t_int i = 0; i < n_pvt_region; ++i)
      {
        boost::array <t_float, 2> dbuf;
        if (reader->read_array (keyword, dbuf) != 2)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        p_ref[i] = dbuf[0];
        rock[i] = dbuf[1];
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  pvt_keywords::pvt_keywords (bs_type_ctor_param)
  {
  }

  pvt_keywords::pvt_keywords (const pvt_keywords &src)
  : bs_refcounter (src), keyword_info_base (src)
  {
  }

  void
  pvt_keywords::register_keywords (sp_objbase &km, std::string provider) const
  {
    BS_SP (keyword_manager_iface) keyword_manager (km, bs_dynamic_cast ());
    BS_ASSERT (keyword_manager);

    keyword_manager->register_keyword ("DENSITY", keyword_handler (DENSITY, 0));
    keyword_manager->register_keyword ("PVTO", keyword_handler (PVTO, 0));
    keyword_manager->register_keyword ("PVDO", keyword_handler (PVDO, 0));
    keyword_manager->register_keyword ("PVTW", keyword_handler (PVTW, 0));
    keyword_manager->register_keyword ("PVDG", keyword_handler (PVDG, 0));
    keyword_manager->register_keyword ("ROCK", keyword_handler (ROCK, 0));
  }


  BLUE_SKY_TYPE_STD_CREATE (pvt_keywords);
  BLUE_SKY_TYPE_STD_COPY (pvt_keywords);
  BLUE_SKY_TYPE_IMPL (pvt_keywords, keyword_info_base, "PVT Keywords", "pvt_keywords", "pvt_keywords");
}
