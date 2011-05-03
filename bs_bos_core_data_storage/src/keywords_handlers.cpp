/**
* @file keywords.h
* @brief functions - keyword handlers
* @author Morozov Andrey
* @date 2008-07-02
*/
#include "bs_bos_core_data_storage_stdafx.h"

//#include "event_base.h"
//#include "event_manager.h"
//#include "fi_params.h"
#include "keyword_manager.h"
#include "localization.h"
#include "constants.h"
#include "data_class.h"


#include BS_FORCE_PLUGIN_IMPORT ()
/*
#include "scal_3p.h"
#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"
*/
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  //BS_TYPE_IMPL_T_MEM(bos_array, amap_strategy_fi::item_array_t)
  //BS_TYPE_IMPL_T_MEM(bos_array, amap_strategy_ii::item_array_t)

//! macro for defining of common variables
#define KH_READER_DEF                                       \
  sp_reader_t reader = params.hdm->get_reader ();           \
  size_t len;                                               \
  len = 0;

  
  keyword_manager::~keyword_manager ()
  {

  }

  
  void keyword_manager::int_array_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_this_t km = params.hdm->get_keyword_manager ();
    sp_pool_t pool = params.hdm->get_pool ();
    
    spv_int this_arr = BS_KERNEL.create_object (v_int::bs_type ());
    t_long ndim = pool->calc_data_dims (keyword);

    this_arr->resize (ndim);
    
    if ((len = reader->read_array (keyword, *this_arr)) != (size_t)ndim)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s") % reader->get_prefix () % keyword);
      }
    pool->set_i_data (keyword, this_arr);
    BOSOUT (section::read_data, level::medium) << "int pool keyword: " << keyword << bs_end;
    BOSOUT (section::read_data, level::medium) << "ndim = " << ndim << bs_end;
  }

  
  void keyword_manager::float_array_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_this_t km = params.hdm->get_keyword_manager ();
    sp_pool_t pool = params.hdm->get_pool ();
    
    //!TODO: check sp!
    spv_float this_arr = BS_KERNEL.create_object (v_float::bs_type ());
    t_long ndim = pool->calc_data_dims (keyword);
    
    this_arr->resize (ndim);

    if ((len = reader->read_array (keyword.c_str(), *this_arr)) != (size_t)ndim)
      {
        bs_throw_exception ((boost::format ("Error in %s: not enough valid arguments for keyword %s") % reader->get_prefix() % keyword).str ());
      }
    pool->set_fp_data (keyword, this_arr);

    BOSOUT (section::read_data, level::medium) << "fp pool keyword: " << keyword << bs_end;
    BOSOUT (section::read_data, level::medium) << "ndim = " << ndim << bs_end;
  }

  void keyword_manager::prop_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_this_t km = params.hdm->get_keyword_manager ();
    sp_idata_t idata = params.hdm->get_data ();
    char buf[CHAR_BUF_LEN] = {0};
    char s_prop[CHAR_BUF_LEN];
    char *start, *end;
    
    t_int i, n, i_prop;
    t_float f_prop;
    std::string &format = km->handlers[keyword].prop_format;
    prop_names_t &names = km->handlers[keyword].prop_names;
    
    n = format.length();
    reader->read_line (buf, CHAR_BUF_LEN);
    start = buf;
    for (i = 0; i <n; ++i)
      {
        if (format[i] == 'i')
        {
          scanf_u (start, &end, &i_prop);
          idata->props->add_property_i (i_prop, names[i], names[i]);
        }
        else if (format[i] == 'f')
        {
          scanf_d (start, &end, &f_prop);
          idata->props->add_property_f (f_prop, names[i], names[i]);
        }
        else if (format[i] == 's')
        {
          scanf_s (start, &end, s_prop);
          idata->props->add_property_s (s_prop, names[i], names[i]);
        }
        else if (format[i] == 'b')
        {
          scanf_s (start, &end, s_prop);
          if (strcmp (s_prop, "YES") == 0 || strcmp (s_prop, "TRUE") == 0 || strcmp (s_prop, "1") == 0 )
            idata->props->add_property_b (1, names[i], names[i]);
          else if (strcmp (s_prop, "NO") == 0 || strcmp (s_prop, "FALSE") == 0 || strcmp (s_prop, "0") == 0 )
            idata->props->add_property_b (0, names[i], names[i]);
          else bs_throw_exception ((boost::format ("Error in %s: not enough valid arguments for keyword %s") % reader->get_prefix() % keyword).str ());
        }
        else
           bs_throw_exception ((boost::format ("Error in %s: not enough valid arguments for keyword %s") % reader->get_prefix() % keyword).str ());

        start = end;
      }


    BOSOUT (section::read_data, level::medium) << "prop keyword: " << keyword << "(" << format << ")" << bs_end;
  }


  /*
  
  void keyword_manager::event_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    typedef smart_ptr <event_manager , true>	sp_event_manager_t;
    sp_event_manager_t em (params.em, bs_dynamic_cast ());
    sp_this_t km = params.hdm->get_keyword_manager ();
    
    
    for (;;)
      {
        reader->read_line (buf, CHAR_BUF_LEN);
        if (sscanf (buf, "%s", key1) != 1)
          {
            bs_throw_exception (boost::format ("Error in %s: bad string (%s)")
              % reader->get_prefix() % buf);
          }
        if (key1[0] == '/')
          {
            em->end_event ();
            break;
          }

        em->process_event (km->current_date, keyword, std::string (buf));
      }
  }
*/

  
  void keyword_manager::TITLE_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    sp_idata_t idata = params.hdm->get_data ();
    
    reader->read_line (buf, CHAR_BUF_LEN);
    idata->props->set_s("title", std::string(buf)); // TODO: , len)
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::OIL_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata = params.hdm->get_data ();

    idata->props->set_b("oil_phase", 1);
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::WATER_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata = params.hdm->get_data ();

    idata->props->set_b("water_phase", 1);
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::GAS_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata = params.hdm->get_data ();
    
    idata->props->set_b("water_phase", 1);
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

/*
  
  void keyword_manager::PROCESS_PARAMS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    
    try
      {
        smart_ptr <fi_params, true> fi_params (params.fi_params, bs_dynamic_cast ());
        fi_params->read_from_keyword_data_file (reader);
      }
    catch (const bs_exception&)
      {
        throw;
      }
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }
*/
  
  void keyword_manager::REPORTS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    char *strt = 0, *end_ptr = 0;
    

    for (;;)
      {
        reader->read_line (buf, CHAR_BUF_LEN);
        if (sscanf (buf, "%s", key1) != 1)
          {
            bs_throw_exception (boost::format ("Error in %s: bad string (%s)")
              % reader->get_prefix() % buf);
          }
        if (key1[0] == '/')
          break;

        strt = buf;
        scanf_s (strt, &end_ptr, key1);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() << ": can't read report type for keyword " << keyword << " from " << strt;
        //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
        //    KH_ASSERT_EXCEPTION
        //  }
        strt = end_ptr;
        int section_level = level::warning;
        if (COMPARE_KEYWORD (key1, "ALL") && COMPARE_KEYWORD (key1, "NO"))
          {
            scanf_u (strt, &end_ptr, &section_level);
          }
        //if (COMPARE_KEYWORD (key1, "ALL")
        //    && COMPARE_KEYWORD (key1, "NO")
        //    && scanf_u (strt, &end_ptr, itmp))
        //  {
        //    out_s << "Error in " << reader->get_prefix() << ": can't read report level for type " << key1 << " for keyword " << keyword << " from " << strt;
        //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
        //    KH_ASSERT_EXCEPTION
        //  }

        // TODO: check section_level value
        if (log::detail::set_section_level_by_name (key1, section_level + level::warning))
          {
            if (std::string (key1) == "ITERS")
              {
                BOSOUT.set_priority (section::solvers, section_level + level::warning);
              }
          }
        else if (!COMPARE_KEYWORD (key1, "ALL"))
          {
            for (int i = section::section_begin; i < section::section_end; ++i)
              {
                BOSOUT.set_priority (i, level::low);
              }
          }
        else if (!COMPARE_KEYWORD (key1, "NO"))
          {
            for (int i = section::section_begin; i < section::section_end; ++i)
              {
                BOSOUT.set_priority (i, level::warning);
              }
          }
        else
          {
            bs_throw_exception (boost::format ("Error in %s: unknown section name %s for keyword %s")
              % reader->get_prefix () % key1 % keyword);
          }
      }
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::REPORTFILE_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    char key2[CHAR_BUF_LEN] = {0};
    

    for (;;)
      {
        reader->read_line (buf, CHAR_BUF_LEN);
        if (buf[0] == '/')
          {
            BOSOUT (section::read_data, level::debug) << "End of REPORT_FILE section" << bs_end;
            break;
          }
        if (sscanf (buf, "%s %s", key1, key2) != 2)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix() % keyword);
          }

        int priority = log::detail::get_log_level_by_name (key2);
        if (priority == -1)
          {
            bs_throw_exception (boost::format ("Error in %s: unknown ('%s') log level (priority) specified for keyword %s")
              % reader->get_prefix () % key2 % keyword);
          }

        if (!log::detail::set_section_level_by_name (key1, priority, BOSOUT.get_channel ()->get_stream ("FILE")))
          {
            bs_throw_exception (boost::format ("Error in %s: unknown ('%s') section name specified for keyword %s")
              %reader->get_prefix () % key1 % keyword);
          }
      }
  }

  
  void keyword_manager::REPORTSCREEN_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    char key2[CHAR_BUF_LEN] = {0};
    

    for (;;)
      {
        reader->read_line (buf, CHAR_BUF_LEN);
        if (buf[0] == '/')
          {
            BOSOUT (section::read_data, level::debug) << "End of REPORT_SCREEN section" << bs_end;
            break;
          }
        if (sscanf (buf, "%s %s", key1, key2) != 2)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix() % keyword);
          }

        int priority = log::detail::get_log_level_by_name (key2);
        if (priority == -1)
          {
            bs_throw_exception (boost::format ("Error in %s: unknown ('%s') log level (priority) specified for keyword %s")
              % reader->get_prefix () % key2 % keyword);
          }

        if (!log::detail::set_section_level_by_name (key1, priority, BOSOUT.get_channel ()->get_stream ("COUT")))
          {
            bs_throw_exception (boost::format ("Error in %s: unknown ('%s') section name specified for keyword %s")
              %reader->get_prefix () % key1 % keyword);
          }
      }
  }

  
  void keyword_manager::STONE1_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata = params.hdm->get_data ();
    
    idata->props->set_i("rpo_model ", STONE1_MODEL);
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::STONE2_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata = params.hdm->get_data ();
    
    idata->props->set_i("rpo_model ", STONE2_MODEL);
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::RELATIVE_PERM_DEFAULT_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata = params.hdm->get_data ();
    
    idata->props->set_i("rpo_model ", RPO_DEFAULT_MODEL);
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  void keyword_manager::SCALECRS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    sp_idata_t idata = params.hdm->get_data ();

    for (;;)
      {
        reader->read_line (buf, CHAR_BUF_LEN);
        if (sscanf (buf, "%s", key1) != 1)
          {
            bs_throw_exception (boost::format ("Error in %s: bad string (%s)")
              % reader->get_prefix() % buf);
          }
        if (key1[0] == '/')
          break;

        if (!COMPARE_KEYWORD (key1, "YES"))
          {
            idata->props->set_b ("scalecrs", 1);
          }
        else if (!COMPARE_KEYWORD (key1, "NO"))  
          {
            idata->props->set_b ("scalecrs", 0);
          }
        else
          {
            bs_throw_exception (boost::format ("Error in %s: unknown argument %s for keyword %s")
              % reader->get_prefix () % key1 % keyword);
          }
        break;  
      }    
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::UNITS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    char key2[CHAR_BUF_LEN] = {0};
    sp_idata_t idata = params.hdm->get_data ();
    
    reader->read_line (buf, CHAR_BUF_LEN);
    // Read and convert UNITS to data
    if (sscanf (buf, "%s", key1) != 1)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }

    int units_in = UNITS_INPUT_DEFAULT;
    if (!COMPARE_KEYWORD (key1, "SI"))
      units_in = UNITS_SI;
    else if (!COMPARE_KEYWORD (key1, "METRIC"))
      units_in = UNITS_METRIC;
    else if (!COMPARE_KEYWORD (key1, "FIELD"))
      units_in = UNITS_FIELD;
    else if (!COMPARE_KEYWORD (key1, "LAB"))
      units_in = UNITS_LAB;
    else
      {
        bs_throw_exception (boost::format ("Error in %s: unknown argument %s for keyword %s")
          % reader->get_prefix () % key1 % keyword);
      }
    // Make output units being equal to input by default
    int units_out = units_in;
    if (sscanf (buf, "%s%s", key1, key2) == 2)
      {
        // two units specified
        if (!COMPARE_KEYWORD (key2, "SI"))
          units_out = UNITS_SI;
        else if (!COMPARE_KEYWORD (key2, "METRIC"))
          units_out = UNITS_METRIC;
        else if (!COMPARE_KEYWORD (key2, "FIELD"))
          units_out = UNITS_FIELD;
        else if (!COMPARE_KEYWORD (key2, "LAB"))
          units_out = UNITS_LAB;
        else
          {
            BOSOUT (section::read_data, level::warning) << "Warning in " << reader->get_prefix()
            << ": trailing garbage " << key2
            << " is ignored for keyword " << keyword << bs_end;
          }
      }

    idata->input_units_converter.set_input_units (units_in);
    idata->input_units_converter.set_output_units (UNITS_INTERNAL_DEFAULT);
    idata->output_units_converter.set_input_units (UNITS_INTERNAL_DEFAULT);
    idata->output_units_converter.set_output_units (units_out);

    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::ROCKCOMP_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    sp_idata_t idata = params.hdm->get_data ();
    t_int rock_region;

    reader->read_line (buf, CHAR_BUF_LEN);
    // Read and convert UNITS to data
    if (sscanf (buf, "%s %d", key1, &rock_region) != 2)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }
    // check read data
    if (rock_region < 1)
      {
        bs_throw_exception (boost::format ("Error in %s: number of rock regions should be greater than 0 for keyword %s")
          % reader->get_prefix () % keyword);
      }
    idata->props->set_i ("rock_region", rock_region);
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::REGNUM_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();

    boost::array <int, 4> regions;
    if ((len = reader->read_array (keyword, regions)) != regions.size ())
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }

    idata->set_region (regions[0], regions[1], regions[2], regions[3]);

    /*
    BS_ASSERT (regions[0] == idata->pvt_region) (regions [0]) (idata->pvt_region);
    BS_ASSERT (regions[1] == idata->sat_region) (regions [1]) (idata->sat_region);
    BS_ASSERT (regions[2] == idata->eql_region) (regions [2]) (idata->eql_region);
    BS_ASSERT (regions[3] == idata->fip_region) (regions [3]) (idata->fip_region);
   */

    BOSOUT (section::read_data, level::medium) <<
    "Keyword " << keyword << ":(" 
      << regions[0] << ", " << regions[1] << ", " 
      << regions[2] << ", " << regions[3] << ")" 
      << bs_end;
  }

  
  void keyword_manager::REGDIMS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    
    boost::array <int, 4> itmp;
    if ((len = reader->read_array (keyword, itmp)) != itmp.size ())
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }

    idata->set_region(itmp[3],itmp[1],itmp[2],itmp[0]);

    BOSOUT (section::read_data, level::medium) <<
    "Keyword " << reader->get_prefix() << ": reading of (" <<
    itmp[0] << ", " << itmp[1] << ", " <<
    itmp[2] << ", " << itmp[3] << ") is successfully" << bs_end;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }


  
  void keyword_manager::EQLDIMS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    
    boost::array <int, 1> itmp;
    if ((len = reader->read_array (keyword, itmp)) != itmp.size ())
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }
    idata->props->set_i("eql_region", itmp[0]);
    if (itmp[0] <= 0)
      {
        bs_throw_exception (boost::format ("Error in %s: number of equilibrium regions in %s must be positive")
          % reader->get_prefix () % keyword);
      }

    idata->equil->resize(EQUIL_TOTAL * itmp[0]); //!TODO:EQUIL_TOTAL instead of 3

    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::TABDIMS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char *strt = 0, *end_ptr = 0;
    sp_idata_t idata = params.hdm->get_data ();
    boost::array <int, 4> itmp;

    reader->read_line (buf, CHAR_BUF_LEN);
    // read number of saturation tables
    strt = buf;
    scanf_u (strt, &end_ptr, &itmp[1]);
    //if (reader->)
    //  {
    //    out_s << "Error in " << reader->get_prefix() << ": can't read number of saturation tables from " << strt;
    //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
    //    KH_ASSERT_EXCEPTION
    //  }

    strt=end_ptr;
    // read number of PVT tables
    scanf_u (strt, &end_ptr, &itmp[0]);
    //if (reader->)
    //  {
    //    out_s << "Error in " << reader->get_prefix() << ": can't read number of PVT tables from " << strt;
    //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
    //    KH_ASSERT_EXCEPTION
    //  }
    strt=end_ptr;
    //skip two currently unesed parameters
    skip_u (strt, &end_ptr, 2);
    //if (reader->)
    //  {
    //    out_s << "Error in " << reader->get_prefix() << ": can't skip currently used parameters from " << strt;
    //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
    //    KH_ASSERT_EXCEPTION
    //  }
    strt = end_ptr;
    // read number of FIP regions
    scanf_u (strt, &end_ptr, &itmp[3]);
    //if (reader->)
    //  {
    //    out_s << "Error in " << reader->get_prefix() << ": can't read the maximum number of FIP regions from " << strt;
    //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
    //    KH_ASSERT_EXCEPTION
    //  }

    // Allocate memory for TABDIMS
    itmp[2] = idata->props->get_i ("eql_region");

    idata->set_region(itmp[0],itmp[1],itmp[2],itmp[3]);

    BOSOUT (section::read_data, level::medium) <<
    "Keyword " << reader->get_prefix() << ": reading of (" <<
    itmp[1] << ", " << itmp[0] << ", 2*, " <<
    itmp[3] << ") is successfull" << bs_end;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  
  void keyword_manager::DENSITY_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    spv_float density = BS_KERNEL.create_object (v_float::bs_type ());
    t_long n_pvt_region = idata->props->get_i("pvt_region");
    
    density->resize (n_pvt_region * 3);
    
    for (t_int i = 0; i < n_pvt_region; i++)
      {
        if ((len = reader->read_array (keyword, *density, 3*i, 3)) != 3)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
      }

    try
      {
        idata -> set_density(density);
      }
    catch (const bs_exception& e)
      {
        bs_throw_exception (boost::format ("Error in %s: %s for keyword %s")
          % reader->get_prefix () % e.what () % keyword);
      }

    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::ROCKTAB_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_rock_region = idata->props->get_i ("rock_region");
    
    if (n_rock_region < 1)
      {
        bs_throw_exception (boost::format ("Error in %s: keyword ROCKCOMP should be used before keyword %s")
          % reader->get_prefix () % keyword);
      }

    idata->rocktab.resize (n_rock_region);

    // Read table for each of region
    std::vector<t_float> dbuf;
    for (t_int i = 0; i < n_rock_region; i++)
      {
        std::vector <t_float> &p_col = idata->rocktab[i].get_column (0);
        std::vector <t_float> &pvm_col = idata->rocktab[i].get_column (1);
        std::vector <t_float> &tm_col = idata->rocktab[i].get_column (2);

        if ((len = reader->read_table (keyword, dbuf, 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        idata->rocktab[i].set_num_rows (t_int (len));
        for (size_t j = 0; j < len; ++j)
          {
            p_col[j]      = dbuf[j * 3 + 0];
            pvm_col[j]    = dbuf[j * 3 + 1];
            tm_col[j]     = dbuf[j * 3 + 2];
          }
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::PVTO_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    int i;
    int j;
    int lj = 0;
    t_float dbuf[DOUB_BUF_LEN] = {0};
    t_float *main_data;
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_pvt_region = idata->props->get_i("pvt_region");
    
    if (!idata->pvto.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for oil has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }

    // Read table for each of region
    for (i = 0; i < n_pvt_region; i++)
      {
        lj = 0;
        char buf[CHAR_BUF_LEN] = {0};
        if (reader->read_line (buf, CHAR_BUF_LEN) <= 0)
          {
            bs_throw_exception (boost::format ("Error in %s: can't read argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        if (sscanf (buf, "%lf%lf%lf%lf", &dbuf[lj * 4],
                    &dbuf[lj * 4 + 1], &dbuf[lj * 4 + 2],
                    &dbuf[lj * 4 + 3]) != 4)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        for (lj = 1;; ++lj)
          {
            if (reader->read_line (buf, CHAR_BUF_LEN) <= 0)
              {
                bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
                  % reader->get_prefix () % keyword);
              }
            if (buf[0] == '/')
              {
                idata->pvto[i].main_data_->resize(lj*4);
                break;
              }

            if ((len = sscanf (buf, "%lf%lf%lf%lf", &dbuf[lj * 4],
                               &dbuf[lj * 4 + 1], &dbuf[lj * 4 + 2],
                               &dbuf[lj * 4 + 3])) != 4)
              {
                if (len == 3)
                  {
                    // 3 values table
                    dbuf[lj * 4 + 3] = dbuf[lj * 4 + 2];
                    dbuf[lj * 4 + 2] = dbuf[lj * 4 + 1];
                    dbuf[lj * 4 + 1] = dbuf[lj * 4];
                    dbuf[lj * 4] = -1.0; //dbuf[(lj - 1) * 4];    // previous
                  }
                else
                  {
                    bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
                      % reader->get_prefix () % keyword);
                  }
              }
          }

        main_data = idata->pvto[i].main_data_->data ();
        // Rows infill
        for (j = 0; j < lj*4; j++)
          {
            main_data[j] = dbuf[j];
          }

        BOSOUT (section::read_data, level::medium) << "lj=" << lj << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::PVDO_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_pvt_region = idata->props->get_i("pvt_region");
    
    if (!idata->pvtdo.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for dead oil has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }

    // Read table for each of region
    for (t_int i = 0; i < n_pvt_region; i++)
      {
        if ((len = reader->read_table (keyword, (*idata->pvtdo[i].main_data_), 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::PVTW_handler(const std::string &keyword, keyword_params_t &params)/********************************/
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_pvt_region = idata->props->get_i("pvt_region");
    
    if (!idata->pvtw.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for water has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }

    // Read table for each of region
    for (t_int i = 0; i < n_pvt_region; i++)
      {
        //idata->pvtw[i].main_data_.resize(4);
        if ((len = reader->read_table (keyword, (*idata->pvtw[i].main_data_), 4)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        if (len != 1)
          {
            bs_throw_exception (boost::format ("Error in %s: you can specify only 4 arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }
      }

    BOSOUT (section::read_data, level::medium) << "pvt_region=" << n_pvt_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::PVDG_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_pvt_region = idata->props->get_i("pvt_region");

    if (!idata->pvtg.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for gas has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }

    // Read table for each of region
    for (t_int i = 0; i < n_pvt_region; i++)
      {
        //idata->pvtg[i].main_data_.resize(DOUB_BUF_LEN);
        if ((len = reader->read_table (keyword, (*idata->pvtg[i].main_data_), 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }


        // Rows infill
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::ROCK_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    boost::array <t_float, 2> dbuf;
    t_float *p_ref = idata->p_ref->data ();
    t_float*rock = idata->rock->data ();
    t_long n_pvt_region = idata->props->get_i("pvt_region");
    
    // Compressibility of rock
    for (t_int i = 0; i < n_pvt_region; ++i)
      {
        if ((len = reader->read_array (keyword, dbuf)) != 2)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        p_ref[i] = dbuf[0];
        rock[i] = dbuf[1];
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::SWOF_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->swof.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SWOF table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
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

  
  void keyword_manager::SGOF_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sgof.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SGOF table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
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


  
  void keyword_manager::SWFN_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->swfn.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SWFN table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
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


  void keyword_manager::SGFN_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sgfn.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SGFN table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
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


  
  void keyword_manager::SOF3_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sof3.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SOF3 table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
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
  
  
  void keyword_manager::SOF2_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    t_long n_scal_region = idata->props->get_i("scal_region");

    if (!idata->sof2.size())
      {
        bs_throw_exception (boost::format ("Error in %s: SOF3 table has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }
  
    // Read table for each of region
    for (t_long i = 0; i < n_scal_region; ++i)
      {
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
  
/*

  
  void keyword_manager::EQUIL_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    boost::array <double, DOUB_BUF_LEN> dbuf;
    char *end_ptr = 0;
    sp_idata_t idata = params.hdm->get_data ();
    sp_this_t km = params.hdm->get_keyword_manager ();

    //! EQUIL keyword enumerated params

    if (!idata->pvto.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for oil has not been initialized yet (keyword: %s")
          % reader->get_prefix () % keyword);
      }
    if (!idata->pvtw.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for water has not been initialized yet (keyword: %s")
          % reader->get_prefix () % keyword);
      }
    if (!idata->pvtg.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for gas has not been initialized yet (keyword: %s")
          % reader->get_prefix () % keyword);
      }
    if (!idata->rock->size())
      {
        bs_throw_exception (boost::format ("Error in %s: rock properties have not been initialized yet (keyword: %s")
          % reader->get_prefix () % keyword);
      }
    if (!idata->equil->size())
      {
        bs_throw_exception (boost::format ("Error in %s: EQUIL table has not been initialized yet (keyword: %s")
          % reader->get_prefix () % keyword);
      }

    if (idata->eql_region == 0)
      idata->eql_region = 1;

    for (size_t i = 0; i < idata->eql_region; ++i)
      {
        reader->read_line (buf, CHAR_BUF_LEN);
        if (buf[0] == '/')
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }

        // Values by default
        dbuf[EQUIL_DAT_DEPTH] = -1.0; //Ddat
        dbuf[EQUIL_DAT_PRESS] = 0.0;  //Pdat
        dbuf[EQUIL_WOC_DEPTH] = -1.0; //Dwoc
        dbuf[EQUIL_WOC_PRESS] = 0.0;  //Pwoc
        dbuf[EQUIL_GOC_DEPTH] = 0.0;  // Dgoc
        dbuf[EQUIL_GOC_PRESS] = 0.0;  // Pgoc
        dbuf[EQUIL_RS_TYPE]   = 0;    // RS type initialization
        dbuf[EQUIL_RV_TYPE]   = 0;    // RV type initialization
        dbuf[EQUIL_NUM_SEC]   = 100;  // number of points
        dbuf[EQUIL_COMP_TYPE] = 0;    // composition type initialization
        dbuf[EQUIL_COMP_ARG]  = 0;    // composition argument initialization

        unwrap (buf);
        scanf_d (buf,     &end_ptr, &dbuf[EQUIL_DAT_DEPTH], 0, 0, "Can't read EQUIL_DEPTH");
        scanf_d (end_ptr, &end_ptr, &dbuf[EQUIL_DAT_PRESS], 0, 0, "Can't read EQUIL_DEPTH_PRESSURE");
        scanf_d (end_ptr, &end_ptr, &dbuf[EQUIL_WOC_DEPTH], 0, 0, "Can't read EQUIL_WOC_DEPTH");
        scanf_d (end_ptr, &end_ptr, &dbuf[EQUIL_WOC_PRESS], 0, 0, "Can't read EQUIL_WOC_DEPTH_PRESSURE");
        scanf_d (end_ptr, &end_ptr, &dbuf[EQUIL_GOC_DEPTH], 0, 0, "Can't read EQUIL_GOC_DEPTH");
        scanf_d (end_ptr, &end_ptr, &dbuf[EQUIL_GOC_PRESS], 0, 0, "Can't read EQUIL_GOC_DEPTH_PRESSURE");
        scanf_d (end_ptr, &end_ptr, &dbuf[EQUIL_RS_TYPE],   0, 0, "Can't read RS initialization type at GOC depth");

        if (dbuf[EQUIL_DAT_DEPTH] < 0) // Ddat
          {
            bs_throw_exception (boost::format ("Error in %s: incorrect value of EQUIL_DEPTH %f")
              % reader->get_prefix () % dbuf[EQUIL_DAT_DEPTH]);
          }

        if (dbuf[EQUIL_DAT_PRESS] < 0) //Pdat
          {
            bs_throw_exception (boost::format ("Error in %s: incorrect pressure value at datum depth %f")
              % reader->get_prefix () % dbuf[EQUIL_DAT_PRESS]);
          }

        if (dbuf[EQUIL_WOC_DEPTH] < 0)
          {
            bs_throw_exception (boost::format ("Error in %s: incorrect value of WOC depth %f")
              % reader->get_prefix () % dbuf[EQUIL_WOC_DEPTH]);
          }

        if (dbuf[EQUIL_GOC_DEPTH] < 0)
          {
            //                rep->print (LOG_READ_SECTION, LOG_ERR,
            //                           GET_TEXT ("Error in %s: incorrect value of GOC depth %lf\n"), r.get_prefix (), dbuf[4]);
            //                  return -8;
            dbuf[EQUIL_GOC_DEPTH] = 0;
          }

        for (size_t ii = 0; ii < EQUIL_TOTAL; ++ii)
          idata->equil[EQUIL_TOTAL * i + ii] = dbuf[ii];
      }

    if (idata->prvd.size())
      idata->prvd.resize(0);
    if (idata->rsvd.size())
      idata->rsvd.resize(0);
    if (idata->pbvd.size())
      idata->pbvd.resize(0);

    idata->init_section = 1;       // flag indicating init section

    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::PRVD_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    
    if (idata->eql_region == 0)
      idata->eql_region = 1;

    // allocate memory for prvd
    idata->prvd.resize(idata->eql_region);

    // Read table for each of region
    shared_vector <double> dbuf;
    for (size_t i = 0; i < idata->eql_region; ++i)
      {
        if ((len =reader->read_table (keyword, dbuf, 2)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }

        idata->prvd[i].set_table_len (i_type_t (len));

        std::vector <double> &dpt = idata->prvd[i].tdepth();
        std::vector <double> &prs = idata->prvd[i].tvalues();
        // Rows infill
        for (size_t j = 0; j < len; ++j)
          {
            dpt[j] = dbuf[j * 2];
            prs[j] = dbuf[j * 2 + 1];
          }

        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::RSVD_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    
    if (idata->eql_region == 0)
      idata->eql_region = 1;

    // allocate memory for rsvd
    idata->rsvd.resize(idata->eql_region);

    // Read table for each of region
    shared_vector <double> dbuf;
    for (size_t i = 0; i < idata->eql_region; ++i)
      {
        if ((len =reader->read_table (keyword, dbuf, 2)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }

        idata->rsvd[i].set_table_len (i_type_t (len));

        std::vector <double> &dpt = idata->rsvd[i].tdepth();
        std::vector <double> &prs = idata->rsvd[i].tvalues();

        dpt.assign (len, 0);
        prs.assign (len, 0);
        for (size_t j = 0; j < len; ++j)
          {
            dpt[j] = dbuf[j * 2 + 0];
            prs[j] = dbuf[j * 2 + 1];
          }

        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  
  void keyword_manager::PBVD_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata = params.hdm->get_data ();
    
    if (idata->eql_region == 0)
      idata->eql_region = 1;

    // allocate memory for rsvd
    idata->pbvd.resize(idata->eql_region);

    // Read table for each of region
    shared_vector <double> dbuf;
    for (size_t i = 0; i < idata->eql_region; ++i)
      {
        if ((len =reader->read_table (keyword, dbuf, 2)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }

        idata->pbvd[i].set_table_len (i_type_t (len));

        std::vector <double> &dpt = idata->pbvd[i].tdepth();
        std::vector <double> &prs = idata->pbvd[i].tvalues();

        dpt.assign (len, 0);
        prs.assign (len, 0);
        for (size_t j = 0; j < len; ++j)
          {
            dpt[j] = dbuf[j * 2 + 0];
            prs[j] = dbuf[j * 2 + 1];
          }

        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  
  void keyword_manager::START_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    sp_this_t km = params.hdm->get_keyword_manager ();
    
    if ((reader->read_line (buf, CHAR_BUF_LEN)) <= 0)
      {
        bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
          % reader->get_prefix () % keyword);
      }
    boost::gregorian::date start=reader->read_date(std::string(buf));
    km->starting_date = boost::posix_time::ptime(start);        // set starting date

    // Current date
    km->current_date = boost::posix_time::ptime(start);;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  
  void keyword_manager::DATE_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    typedef event_manager  event_manager_t;
    typedef smart_ptr <event_manager_t, true>	sp_event_manager_t;
    sp_event_manager_t em (params.em, bs_dynamic_cast ());
    sp_this_t km = params.hdm->get_keyword_manager ();
    
    if ((reader->read_line (buf, CHAR_BUF_LEN)) <= 0)
      {
        bs_throw_exception (boost::format ("Error in %s: can't read arguments for keyword %s")
          % reader->get_prefix () % keyword);
      }

    km->current_date=boost::posix_time::ptime(reader->read_date (std::string(buf)));
    if (em->event_list.find (km->current_date)==em->event_list.end ())
      {
        std::list <typename event_manager_t::sp_event_base> tl;
        em->event_list.insert (std::make_pair (km->current_date, tl));
      }

    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::DATES_handler(const std::string &keyword, keyword_params_t &params)
  {
    DATE_handler(keyword, params);
  }

  
  void keyword_manager::TSTEP_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    shared_vector <double> dtmp;
    dtmp.resize(MAX_TIME_STEPS_DEF);
    typedef event_manager  event_manager_t;
    typedef smart_ptr <event_manager_t, true>	sp_event_manager_t;
    sp_event_manager_t em (params.em, bs_dynamic_cast ());
    sp_this_t km = params.hdm->get_keyword_manager ();

    if ((len = reader->read_array (keyword, dtmp)) < 1)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix () % keyword);
      }

    //!
    //! dtmp[0] <- days
    //! dtmp[0]*24*60*60 <-seconds
    //! dtmp[0]*24*60*60*1000 <-millisecs

    for (size_t i = 0; i < (size_t) len; ++i)
      {
        km->current_date += boost::posix_time::millisec (dtmp[i]*24*60*60*1000);

        if (em->event_list.find (km->current_date) == em->event_list.end ())
          {
            std::list <typename event_manager_t::sp_event_base> tl;
            em->event_list.insert (std::make_pair (km->current_date, tl));
          }
      }

    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  
  void keyword_manager::TSTEPS_handler(const std::string &keyword, keyword_params_t &params)
  {
    TSTEP_handler(keyword, params);
  }
  */
//bs stuff

  /*
  KH_SPEC(int_array_handler)
  KH_SPEC(float_array_handler)
//  KH_SPEC(event_handler)
  KH_SPEC(TITLE_handler)
  KH_SPEC(OIL_handler)
  KH_SPEC(WATER_handler)
  KH_SPEC(GAS_handler)
  //KH_SPEC(PROCESS_PARAMS_handler)
  KH_SPEC(RESTART_handler)
  KH_SPEC(REPORTS_handler)
  KH_SPEC(REPORTFILE_handler)
  KH_SPEC(REPORTSCREEN_handler)
  KH_SPEC(STONE1_handler)
  KH_SPEC(STONE2_handler)
  KH_SPEC(RELATIVE_PERM_DEFAULT_handler)
  KH_SPEC(SCALECRS_handler)
  KH_SPEC(UNITS_handler)
  //KH_SPEC(DIMENS_handler)
  KH_SPEC(ROCKCOMP_handler)
  KH_SPEC(REGDIMS_handler)
  KH_SPEC(REGNUM_handler)
  KH_SPEC(EQLDIMS_handler)
  KH_SPEC(TABDIMS_handler)
  //KH_SPEC(COORD_handler)
  //KH_SPEC(ZCORN_handler)
  KH_SPEC(MINPV_handler)
  KH_SPEC(MINSV_handler)
  KH_SPEC(DENSITY_handler)
  KH_SPEC(ROCKTAB_handler)
  KH_SPEC(PVTO_handler)
  KH_SPEC(PVDO_handler)
  KH_SPEC(PVTW_handler)
  KH_SPEC(PVDG_handler)
  KH_SPEC(ROCK_handler)
  
  KH_SPEC(SWOF_handler)
  KH_SPEC(SGOF_handler)
  KH_SPEC(SWFN_handler)
  KH_SPEC(SGFN_handler)
  KH_SPEC(SOF3_handler)
  KH_SPEC(SOF2_handler)
  KH_SPEC(EQUIL_handler)
  KH_SPEC(PRVD_handler)
  KH_SPEC(RSVD_handler)
  KH_SPEC(PBVD_handler)
  KH_SPEC(START_handler)
  KH_SPEC(DATE_handler)
  KH_SPEC(DATES_handler)
  KH_SPEC(TSTEP_handler)
  KH_SPEC(TSTEPS_handler)
*/

//  KH_SPEC(WELLDIMS_handler)

  //!TODO: kill next string after debug
  //template class keyword_manager <base_strategy_di>;
  //template class keyword_manager <base_strategy_fi>;

}//nsbs
