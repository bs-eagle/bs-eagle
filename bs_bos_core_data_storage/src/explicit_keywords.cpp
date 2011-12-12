/**
 *       \file  explicit_keywords.cpp
 *      \brief  Keywords for EXPLICIT model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  29.04.2011
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#include "bs_bos_core_data_storage_stdafx.h"
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif
#include "explicit_keywords.hpp"
#include "keyword_manager_iface.h"
#include "init_model_iface.hpp"
#include "read_class.h"
#include "data_class.h"

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
    PRVD (std::string const &keyword, keyword_params &params)
    {
      BS_SP (FRead) reader = params.hdm->get_reader ();
      BS_SP (idata) idata = params.hdm->get_data ();

      t_long eql_region = idata->props->get_i ("eql_region");
      if (eql_region == 0)
        {
          bs_throw_exception (boost::format ("Error in %s: eql_region == 0 (keyword: %s)")
            % reader->get_prefix () % keyword);
        }

      idata->prvd.resize(eql_region);

      // Read table for each of region
      std::vector <double> dbuf;
      for (t_long i = 0; i < eql_region; ++i)
        {
          t_long len = reader->read_table (keyword, dbuf, 2);
          if (len < 1)
            {
              bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
                % reader->get_prefix () % keyword);
            }

          idata->prvd[i].set_table_len (len);

          std::vector <double> &dpt = idata->prvd[i].tdepth();
          std::vector <double> &prs = idata->prvd[i].tvalues();
          // Rows infill
          for (t_long j = 0; j < len; ++j)
            {
              dpt[j] = dbuf[j * 2];
              prs[j] = dbuf[j * 2 + 1];
            }

          BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
        }
      BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
    }

    void
    activate (std::string const &, keyword_params &params)
    {
      BS_SP (keyword_manager_iface) keyword_manager = params.hdm->get_keyword_manager ();
      BS_ASSERT (keyword_manager);

      BS_SP (init_model_iface) init_model (BS_KERNEL.create_object ("explicit_init_model"), bs_dynamic_cast ());
      params.hdm->set_init_model (init_model);
      params.hdm->get_prop()->add_property_i(0, "init", "init type");
      // FIXME: npy_intp
      npy_intp dimens[] = {1, 0, 1, 0, 1, 0};
      int n_phases;

      n_phases = params.hdm->get_prop()->get_b("oil_phase");
      n_phases += params.hdm->get_prop()->get_b("water_phase");
      n_phases += params.hdm->get_prop()->get_b("gas_phase");

      keyword_manager->register_fp_pool_keyword ("PRESSURE", dimens, 200.0, 0);

      if (n_phases > 1)
        {
          if (params.hdm->get_prop()->get_b("water_phase"))
            keyword_manager->register_fp_pool_keyword ("SWAT",    dimens, 0.3, 0);

          if (params.hdm->get_prop()->get_b("gas_phase"))
            {
              keyword_manager->register_fp_pool_keyword ("SGAS",    dimens, 0, 0);
              keyword_manager->register_fp_pool_keyword ("SOIL",    dimens, 0, 0);
              keyword_manager->register_fp_pool_keyword ("RS",      dimens, 0, 0);
              keyword_manager->register_fp_pool_keyword ("PBUB",    dimens, 0, 0);
            }
        }
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
