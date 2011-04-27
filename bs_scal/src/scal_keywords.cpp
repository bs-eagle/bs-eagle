/**
 *       \file  scal_keywords.cpp
 *      \brief  keywords for SCAL
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  27.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "bs_scal_stdafx.h"
#include "scal_keywords.hpp"
#include "keyword_manager_iface.h"
#include "data_class.h"
#include "read_class.h"
#include "scal_3p_iface.hpp"

namespace blue_sky 
{
  typedef void (*handler_callback) (BS_SP (scal_3p_iface), stdv_float const &, t_long);

  void
  handler (std::string const &keyword, keyword_params &params, handler_callback callback)
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    BS_SP (scal_3p_iface) scal = params.hdm->get_scal ();

    t_long regions = idata->props->get_i ("sat_region");
    stdv_float data;
    for (t_long i = 0; i < regions; ++i)
      {
        if (reader->read_table (keyword, data, 4) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough arguments for keyword %s")
                                % reader->get_prefix () % keyword);
          }

        callback (scal, data, i);
        data.clear ();
      }
  }

  void
  SWOF_callback (BS_SP (scal_3p_iface) scal, stdv_float const &data, t_long)
  {
    spv_float d = BS_KERNEL.create_object (v_float::bs_type ());
    d->init (data);
    scal->get_water_data ()->add_spof (d, true);
  }
  void
  SGOF_callback (BS_SP (scal_3p_iface) scal, stdv_float const &data, t_long)
  {
    spv_float d = BS_KERNEL.create_object (v_float::bs_type ());
    d->init (data);
    scal->get_gas_data ()->add_spof (d, false);
  }

  void
  SWFN_callback (BS_SP (scal_3p_iface) scal, stdv_float const &data, t_long i)
  {
    spv_float d = BS_KERNEL.create_object (v_float::bs_type ());
    d->init (data);
    scal->get_water_data ()->add_spfn (d, i, true);
  }
  void
  SGFN_callback (BS_SP (scal_3p_iface) scal, stdv_float const &data, t_long i)
  {
    spv_float d = BS_KERNEL.create_object (v_float::bs_type ());
    d->init (data);
    scal->get_gas_data ()->add_spfn (d, i, false);
  }

  void
  SOF3_callback (BS_SP (scal_3p_iface) scal, stdv_float const &data, t_long i)
  {
    spv_float d = BS_KERNEL.create_object (v_float::bs_type ());
    d->init (data);
    scal->get_water_data ()->add_sof3 (d, i, true);
    scal->get_gas_data ()->add_sof3 (d, i, false);
  }

  void
  SWOF (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SWOF_callback);
  }
  void
  SGOF (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SGOF_callback);
  }
  void
  SWFN (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SWFN_callback);
  }
  void
  SGFN (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SGFN_callback);
  }
  void
  SOF3 (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SOF3_callback);
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
    keyword_manager->register_keyword ("SOF3", keyword_handler (SOF3, 0));

    // FIXME: npy_intp
    npy_intp dimens[] = {1, 0, 1, 0, 1, 0};
    keyword_manager->register_fp_pool_keyword ("SGL",   dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("SGU",   dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("SOGCR", dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("SGCR",  dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("SWL",   dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("SWU",   dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("SOWCR", dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("SWCR",  dimens, 0, 0);
  }

  BLUE_SKY_TYPE_STD_CREATE (scal_keywords);
  BLUE_SKY_TYPE_STD_COPY (scal_keywords);
  BLUE_SKY_TYPE_IMPL (scal_keywords, keyword_info_base, "SCAL Keywords", "scal_keywords", "scal_keywords");

}
