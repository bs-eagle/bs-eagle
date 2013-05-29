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
#include "bos_reader_iface.h"
#include "scal_3p_iface.hpp"

namespace blue_sky 
{
  typedef void (*handler_callback) (BS_SP (scal_3p_iface), t_float *, t_int, t_int, t_int);

  void
  handler (std::string const &keyword, keyword_params &params, handler_callback callback, const t_int n_columns)
  {
    BS_SP (bos_reader_iface) reader = params.hdm->get_reader ();
    BS_SP (idata) idata = params.hdm->get_data ();
    BS_SP (scal_3p_iface) scal = params.hdm->get_scal ();

    t_int regions = idata->props->get_i (L"sat_region");
    t_int n_rows;
    
    //spv_float data = BS_KERNEL.create_object (v_float::bs_type ());
    std::vector<t_float> d;
    d.resize (4096);

    for (t_long i = 0; i < regions; ++i)
      {
        if ((n_rows = reader->read_fp_table (keyword.c_str (), &d[0], 4096, n_columns)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough arguments for keyword %s")
                                % reader->get_prefix () % keyword);
          }

        callback (scal, &d[0], i, n_rows, n_columns);
      }
  }

  void
  SWOF_callback (BS_SP (scal_3p_iface) scal, t_float *data, t_int region_index, t_int n_rows, t_int n_cols)
  {
    BS_SP (table_iface) tbl = scal->get_table (FI_PHASE_WATER, region_index);
    tbl->convert_from_tf_array(n_rows, n_cols, data);
  }
  void
  SGOF_callback (BS_SP (scal_3p_iface) scal, t_float *data, t_int region_index, t_int n_rows, t_int n_cols)
  {
    BS_SP (table_iface) tbl = scal->get_table (FI_PHASE_GAS, region_index);
    tbl->convert_from_tf_array(n_rows, n_cols, data);
  }

  void
  SWFN_callback (BS_SP (scal_3p_iface) scal, t_float *data, t_int region_index, t_int n_rows, t_int n_cols)
  {
    BS_SP (table_iface) tbl = scal->get_table (FI_PHASE_WATER, region_index);
    tbl->convert_from_tf_array(n_rows, n_cols, data);
  }
  void
  SGFN_callback (BS_SP (scal_3p_iface) scal, t_float *data, t_int region_index, t_int n_rows, t_int n_cols)
  {
    BS_SP (table_iface) tbl = scal->get_table (FI_PHASE_GAS, region_index);
    tbl->convert_from_tf_array(n_rows, n_cols, data);
  }

  void
  SOF3_callback (BS_SP (scal_3p_iface) scal, t_float *data, t_int region_index, t_int n_rows, t_int n_cols)
  {
    BS_SP (table_iface) tbl = scal->get_table (FI_PHASE_OIL, region_index);
    tbl->convert_from_tf_array(n_rows, n_cols, data);
  }

  void
  SOF2_callback (BS_SP (scal_3p_iface) scal, t_float *data, t_int region_index, t_int n_rows, t_int n_cols)
  {
    // TODO: sof2 for 2-phase models only => need check for phases
    // currently write this table for both water and gas data and use one of them (depends on phase which is present)
    BS_SP (table_iface) tbl = scal->get_table (FI_PHASE_OIL, region_index);
    tbl->convert_from_tf_array(n_rows, n_cols, data);
  }
  
 void
  SWOF (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SWOF_callback, SPOF_KEYWORD_COLUMNS);
  }
  void
  SGOF (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SGOF_callback, SPOF_KEYWORD_COLUMNS);
  }
  void
  SWFN (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SWFN_callback, SPFN_KEYWORD_COLUMNS);
  }
  void
  SGFN (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SGFN_callback, SPFN_KEYWORD_COLUMNS);
  }
  void
  SOF3 (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SOF3_callback, SOF3_KEYWORD_COLUMNS);
  }

  void
  SOF2 (std::string const &keyword, keyword_params &params)
  {
    handler (keyword, params, SOF2_callback, SOF2_KEYWORD_COLUMNS);
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
    keyword_manager->register_keyword ("SOF2", keyword_handler (SOF2, 0));

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
    keyword_manager->register_fp_pool_keyword ("PCW",   dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("PCG",   dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("KRO",   dimens, 0, 0);
    keyword_manager->register_fp_pool_keyword ("KRW",   dimens, 0, 0);
  }

  BLUE_SKY_TYPE_STD_CREATE (scal_keywords);
  BLUE_SKY_TYPE_STD_COPY (scal_keywords);
  BLUE_SKY_TYPE_IMPL (scal_keywords, keyword_info_base, "SCAL Keywords", "scal_keywords", "scal_keywords");

}
