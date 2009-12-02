/**
* @file keywords.h
* @brief functions - keyword handlers
* @author Morozov Andrey
* @date 2008-07-02
*/
#include "stdafx.h"

#include "event_base.h"
#include "event_manager.h"
#include "fi_params.h"
#include "keyword_manager.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "scal_3p.h"
#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"
#include "data_class.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  //BS_TYPE_IMPL_T_MEM(bos_array, amap_strategy_fi::item_array_t)
  //BS_TYPE_IMPL_T_MEM(bos_array, amap_strategy_ii::item_array_t)

//! macro for defining of common variables
#define KH_READER_DEF                                       \
  sp_reader_t reader (params.reader, bs_dynamic_cast ());   \
  size_t len;                                               \
  len = 0;

  template <typename strategy_t>
  keyword_manager<strategy_t>::~keyword_manager ()
  {

  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::int_array_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_this_t km (params.km, bs_dynamic_cast ());
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    index_t* dimens, def_value;
    //!TODO: check sp!
    array_uint8_t this_arr;
    int arr_index = km->handlers[keyword].index_in_pool;
    
    if (arr_index >= 0)
      {
        dimens = &i_pool_sizes[ARRAY_POOL_TOTAL*arr_index];
        def_value = i_pool_default_values[arr_index]; 
        if (!idata->i_map->contain(arr_index))
          {
            idata->i_map->create_item (arr_index, dimens, def_value);
          }
       }
    else
      {
        if (arr_index != -1)
          {
            bs_throw_exception (boost::format ("Wrong handler (%d) invoked for pool keyword %s")
              % arr_index % keyword);
          }
        arr_index = ARR_TOTAL_INT;   
        while (idata->i_map->contain(arr_index))
          {
            arr_index++;
          }
        def_value = km->handlers[keyword].int_def_value;
        dimens = km->handlers[keyword].dimens;
        idata->i_map->create_item (arr_index, dimens, def_value);
        idata->ahelper.add_correspondence_i (keyword, arr_index);
      }
    
    
    this_arr = tools::get_non_empty ((*idata->i_map)[arr_index].array);

    index_t nx = 0, ny = 0, nz = 0, ndim = 0;
    idata->i_map->get_dimens (arr_index, nx, ny, nz, ndim);
    if ((len = reader->read_array (keyword, this_arr)) != (size_t)ndim)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s") % reader->get_prefix () % keyword);
      }
    
    // Check array validity: all elements should be in 1...eql_region
    if (!COMPARE_KEYWORD (keyword.c_str(), "ROCKNUM"))
      for (index_t i = 0; i < ndim; i++)
        {
          if (this_arr[i] < 1 || this_arr[i] > idata->rock_region)
            {
              bs_throw_exception (boost::format ("Error in %s: value %f of element %d of array %s is out of range")
                % reader->get_prefix () % this_arr[i] % i % keyword);
            }
        }
    // Check array validity: all elements should be in 1...eql_region
    if (!COMPARE_KEYWORD (keyword.c_str(), "EQLNUM"))
      for (index_t i = 0; i < ndim; i++)
        {
          if (this_arr[i] < 1 || this_arr[i] > idata->eql_region)
            {
              bs_throw_exception (boost::format ("Error in %s: value %f of element %d of array %s is out of range")
                % reader->get_prefix () % this_arr[i] % i % keyword);
            }
        }
    // Check array validity: all elements should be in 1...sat_region
    if (!COMPARE_KEYWORD (keyword.c_str(), "SATNUM"))
      for (index_t i = 0; i < ndim; i++)
        {
          if (this_arr[i] < 1 || this_arr[i] > idata->sat_region)
            {
              bs_throw_exception (boost::format ("Error in %s: value %f of element %d of array %s is out of range")
                % reader->get_prefix () % this_arr[i] % i % keyword);
            }
        }
    // Check array validity: all elements should be in 1...pvt_region
    if (!COMPARE_KEYWORD (keyword.c_str(),"PVTNUM"))
      for (index_t i = 0; i < ndim; i++)
        {
          if (this_arr[i] < 1 || this_arr[i] > idata->pvt_region)
            {
              bs_throw_exception (boost::format ("Error in %s: value %f of element %d of array %s is out of range")
                % reader->get_prefix () % this_arr[i] % i % keyword);
            }
        }
    // Check array validity: all elements should be 0 or 1
    if (!COMPARE_KEYWORD (keyword.c_str(), "ACTNUM"))
      for (index_t i = 0; i < ndim; i++)
        {
          if (this_arr[i] != 0)
            this_arr[i] = 1;
        }
    // Check array validity: all elements should be in 1...fip_region
    if (!COMPARE_KEYWORD (keyword.c_str(), "FIPNUM"))
      for (index_t i = 0; i < ndim; i++)
        {
          if (this_arr[i] < 1 || this_arr[i] > idata->fip_region)
            {
              bs_throw_exception (boost::format ("Error in %s: value %f of element %d of array %s is out of range")
                % reader->get_prefix () % this_arr[i] % i % keyword);
            }
        }
    // Check array validity: all elements should be 1 first type, 2 second type or 3 read from file
    if (!COMPARE_KEYWORD (keyword.c_str(), "BNDNUM"))
      for (index_t i = 0; i < ndim; i++)
        {
          if (this_arr[i] < 1 && this_arr[i] > 3)
            {
              index_t i_c, j_c, k_c;
              index_t c_n = i;
              k_c = c_n / (idata->dimens.nx * idata->dimens.ny);
              c_n -= k_c * (idata->dimens.nx * idata->dimens.ny);
              j_c = c_n / idata->dimens.nx;
              i_c = c_n - j_c * idata->dimens.nx;

              BOSWARN (section::read_data, level::warning) << "Warning: unknown boundary condition in cell ["
              << i_c << " , " << j_c << " , " << k_c << "]" << bs_end;
              this_arr[i] = 2;
            }
        }
    BOSOUT (section::read_data, level::medium) 
      << "Keyword: " << keyword << "\n"
      << "ndim = " << ndim << bs_end;
    // launch second handler if any
    if (km->handlers[keyword].second_handle_function)
      km->handlers[keyword].second_handle_function (keyword, params);
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::float_array_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_this_t km (params.km, bs_dynamic_cast ());
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    //!TODO: check sp!
    array_float16_t this_arr;
    int arr_index = km->handlers[keyword].index_in_pool;
    index_t* dimens;
    item_t def_value;
    
    if (arr_index >= 0)
      {
        dimens = &d_pool_sizes[ARRAY_POOL_TOTAL*arr_index];
        def_value = d_pool_default_values[arr_index]; 
        BOSOUT (section::read_data, level::debug) <<
          boost::format ("[%s] default values: %f\n\tdimens: [%d, %d, %d, %d, %d, %d]") 
          % keyword % def_value 
          % dimens[0] % dimens[1] %dimens[2]
          % dimens[3] % dimens[4] %dimens[5]
          << bs_end;

        if (!dimens)
          {
            BOSOUT (section::read_data, level::debug) <<
              "dimens is null" << bs_end;
          }

        if (!idata->d_map->contain(arr_index))
          {
            idata->d_map->create_item (arr_index, dimens, def_value);
          }
        else
          {
            BOSOUT (section::read_data, level::debug) <<
              boost::format ("Keyword %s alread in pool") % keyword << bs_end;
          }
       }
    else
      {
        if (arr_index != -2)
          {
            bs_throw_exception (boost::format ("Wrong handler (%d) invoked for pool keyword %s")
              % arr_index % keyword);
          }
        arr_index = ARR_TOTAL_DOUBLE;   
        while (idata->d_map->contain(arr_index))
          {
            arr_index++;
          }
        def_value = km->handlers[keyword].float_def_value;
        dimens = km->handlers[keyword].dimens;
        idata->d_map->create_item (arr_index, dimens, def_value);
        idata->ahelper.add_correspondence_d (keyword, arr_index);
      }

    this_arr = tools::get_non_empty ((*idata->d_map)[arr_index].array);

    index_t nx = 0, ny = 0, nz = 0, ndim = 0;
    idata->d_map->get_dimens (arr_index, nx, ny, nz, ndim);

    if ((len = reader->read_array (keyword.c_str(), this_arr)) != (size_t)ndim)
      {
        bs_throw_exception ((boost::format ("Error in %s: not enough valid arguments for keyword %s") % reader->get_prefix() % keyword).str ());
      }
    BOSOUT (section::read_data, level::medium) << "Keyword: " << keyword << bs_end;
    BOSOUT (section::read_data, level::medium) << "ndim = " << ndim << bs_end;
    
    // launch second handler if any
    if (km->handlers[keyword].second_handle_function)
      km->handlers[keyword].second_handle_function (keyword, params);
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::event_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    typedef smart_ptr <event_manager <strategy_t>, true>	sp_event_manager_t;
    sp_event_manager_t em (params.em, bs_dynamic_cast ());
    sp_this_t km (params.km, bs_dynamic_cast ());
    
    
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::TITLE_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    reader->read_line (buf, CHAR_BUF_LEN);
    idata->title = std::string(buf); // TODO: , len)
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::OIL_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (!(idata->fi_phases & (1 << FI_PHASE_OIL)))
      ++idata->fi_n_phase;
    idata->fi_phases |= 1 << FI_PHASE_OIL;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::WATER_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (!(idata->fi_phases & (1 << FI_PHASE_WATER)))
      ++idata->fi_n_phase;
    idata->fi_phases |= 1 << FI_PHASE_WATER;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::GAS_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (!(idata->fi_phases & (1 << FI_PHASE_GAS)))
      ++idata->fi_n_phase;
    idata->fi_phases |= 1 << FI_PHASE_GAS;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::PROCESS_PARAMS_handler(const std::string &keyword, keyword_params_t &params)
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::REPORTS_handler(const std::string &keyword, keyword_params_t &params)
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::REPORTFILE_handler(const std::string &keyword, keyword_params_t &params)
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::REPORTSCREEN_handler(const std::string &keyword, keyword_params_t &params)
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::STONE1_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    idata->rpo_model = STONE1_MODEL;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::STONE2_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    idata->rpo_model = STONE2_MODEL;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::RELATIVE_PERM_DEFAULT_handler(const std::string &keyword, keyword_params_t &params)
  {
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    idata->rpo_model = RPO_DEFAULT_MODEL;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::UNITS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    char key2[CHAR_BUF_LEN] = {0};
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::ROCKCOMP_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char key1[CHAR_BUF_LEN] = {0};
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    reader->read_line (buf, CHAR_BUF_LEN);
    // Read and convert UNITS to data
    if (sscanf (buf, "%s %d", key1, &idata->rock_region.data ()) != 2)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }
    // check read data
    if (idata->rock_region < 1)
      {
        bs_throw_exception (boost::format ("Error in %s: number of rock regions should be greater than 0 for keyword %s")
          % reader->get_prefix () % keyword);
      }
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::REGNUM_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    shared_vector<int> itmp;
    itmp.resize(4);
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    if ((len = reader->read_array (keyword, itmp)) != 4)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }

    idata->set_region(itmp[0],itmp[1],itmp[2],itmp[3]);

    BOSOUT (section::read_data, level::medium) <<
    "Keyword " << keyword << ":(" <<
    idata->pvt_region << ", " << idata->sat_region << ", " <<
    idata->eql_region << ", " << idata->fip_region << ")" << bs_end;

    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::REGDIMS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    shared_vector<int> itmp;
    itmp.resize(4);
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    

    if ((len = reader->read_array (keyword, itmp)) != 4)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }

    idata->set_region(itmp[3],itmp[1],itmp[2],itmp[0]);

    BOSOUT (section::read_data, level::medium) <<
    "Keyword " << reader->get_prefix() << ": reading of (" <<
    idata->pvt_region << ", " << idata->sat_region << ", " <<
    idata->eql_region << ", " << idata->fip_region << ") is successfully" << bs_end;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }


  template <class strategy_t>
  void keyword_manager<strategy_t>::EQLDIMS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    shared_vector <int> itmp;
    itmp.resize(1);
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    

    if ((len = reader->read_array (keyword, itmp)) != 1)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }
    idata->eql_region=itmp[0];
    if (idata->eql_region <= 0)
      {
        bs_throw_exception (boost::format ("Error in %s: number of equilibrium regions in %s must be positive")
          % reader->get_prefix () % keyword);
      }

    idata->equil.resize(EQUIL_TOTAL * idata->eql_region); //!TODO:EQUIL_TOTAL instead of 3

    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::TABDIMS_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    char *strt = 0, *end_ptr = 0;
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    reader->read_line (buf, CHAR_BUF_LEN);
    // read number of saturation tables
    strt = buf;
    scanf_u (strt, &end_ptr, &idata->sat_region.data ());
    //if (reader->)
    //  {
    //    out_s << "Error in " << reader->get_prefix() << ": can't read number of saturation tables from " << strt;
    //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
    //    KH_ASSERT_EXCEPTION
    //  }

    strt=end_ptr;
    // read number of PVT tables
    scanf_u (strt, &end_ptr, &idata->pvt_region.data ());
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
    scanf_u (strt, &end_ptr, &idata->fip_region.data ());
    //if (reader->)
    //  {
    //    out_s << "Error in " << reader->get_prefix() << ": can't read the maximum number of FIP regions from " << strt;
    //    //BOSERR (section::read_data, level::error) << out_s << bs_end;
    //    KH_ASSERT_EXCEPTION
    //  }

    // Allocate memory for TABDIMS
    if (idata->eql_region < 1)
      idata->eql_region = 1;

    idata->set_region (idata->pvt_region, idata->fip_region, idata->sat_region, idata->eql_region);

    BOSOUT (section::read_data, level::medium) <<
    "Keyword " << reader->get_prefix() << ": reading of (" <<
    idata->pvt_region << ", " << idata->sat_region << ", " <<
    idata->eql_region << ", " << idata->fip_region << ") is successfully" << bs_end;
    BOSOUT (section::read_data, level::medium) << keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::MINPV_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    shared_vector <double> tmp;
    tmp.resize(1);
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    
    if ((len = reader->read_array (keyword, tmp, 0)) != 1)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }
    idata->minimal_pore_volume=tmp[0];
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }
  
  template <class strategy_t>
  void keyword_manager<strategy_t>::MINSV_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    shared_vector <double> tmp;
    tmp.resize(1);
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if ((len = reader->read_array (keyword, tmp, 0)) != 1)
      {
        bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
          % reader->get_prefix() % keyword);
      }
    idata->minimal_splice_volume=tmp[0];
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::DENSITY_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    shared_vector <double> density;
    
    density.resize (idata->pvt_region*3);
    
    for (index_t i = 0; i < (int) idata->pvt_region; i++)
      {
        if ((len = reader->read_array (keyword, density, 3*i, 3)) != 3)
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::ROCKTAB_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (idata->rock_region < 1)
      {
        bs_throw_exception (boost::format ("Error in %s: keyword ROCKCOMP should be used before keyword %s")
          % reader->get_prefix () % keyword);
      }

    idata->rocktab.resize (idata->rock_region);

    // Read table for each of region
    shared_vector <double> dbuf;
    for (index_t i = 0; i < idata->rock_region; i++)
      {
        std::vector <float> &p_col = idata->rocktab[i].get_column (0);
        std::vector <float> &pvm_col = idata->rocktab[i].get_column (1);
        std::vector <float> &tm_col = idata->rocktab[i].get_column (2);

        if ((len = reader->read_table (keyword, dbuf, 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        idata->rocktab[i].set_num_rows (index_t (len));
        for (size_t j = 0; j < len; ++j)
          {
            p_col[j]      = dbuf[j * 3 + 0];
            pvm_col[j]    = dbuf[j * 3 + 1];
            tm_col[j]     = dbuf[j * 3 + 2];
          }
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::PVTO_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    int i;
    int j;
    int lj = 0;
    double dbuf[DOUB_BUF_LEN] = {0};
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (!idata->pvto.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for oil has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }

    // Read table for each of region
    for (i = 0; i < (int) idata->pvt_region; i++)
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
                idata->pvto[i].main_data_.resize(lj*4);
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

        // Rows infill
        for (j = 0; j < lj*4; j++)
          {
            idata->pvto[i].main_data_[j] = dbuf[j];
          }

        BOSOUT (section::read_data, level::medium) << "lj=" << lj << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::PVDO_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (!idata->pvtdo.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for dead oil has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }

    // Read table for each of region
    for (index_t i = 0; i < idata->pvt_region; i++)
      {
        if ((len = reader->read_table (keyword, idata->pvtdo[i].main_data_, 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::PVTW_handler(const std::string &keyword, keyword_params_t &params)/********************************/
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (!idata->pvtw.size())
      {
        bs_throw_exception (boost::format ("Error in %s: PVT table for water has not been initialized yet (keyword: %s)")
          % reader->get_prefix () % keyword);
      }

    // Read table for each of region
    for (index_t i = 0; i < idata->pvt_region; i++)
      {
        //idata->pvtw[i].main_data_.resize(4);
        if ((len = reader->read_table (keyword, idata->pvtw[i].main_data_, 4)) < 1)
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

    BOSOUT (section::read_data, level::medium) << "pvt_region=" << idata->pvt_region << bs_end;
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::PVDG_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    // Read table for each of region
    for (index_t i = 0; i < idata->pvt_region; i++)
      {
        //idata->pvtg[i].main_data_.resize(DOUB_BUF_LEN);
        if ((len = reader->read_table (keyword, idata->pvtg[i].main_data_, 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }


        // Rows infill
        BOSOUT (section::read_data, level::medium) << "len=" << len << " i=" << i << bs_end;
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::ROCK_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    boost::array <item_t, 2> dbuf;

    // Compressibility of rock
    for (index_t i = 0; i < idata->pvt_region; ++i)
      {
        if ((len = reader->read_array (keyword, dbuf)) != 2)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }
        idata->p_ref[i] = dbuf[0];
        idata->rock[i] = dbuf[1];
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  template <class strategy_t>
  void keyword_manager<strategy_t>::SWOF_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    typename strategy_t::item_array_t swof;
    typedef smart_ptr <scal_3p <strategy_t>, true>	sp_scal_3p_t;
    sp_scal_3p_t scal_3p (params.scal_3p);
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    // Read table for each of region
    for (index_t i = 0; i < idata->sat_region; ++i)
      {
        if ((len = reader->read_table (keyword, swof, 4)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        scal_3p->get_water_data()->add_spof (swof, true);
        swof.clear ();
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::SGOF_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    typename strategy_t::item_array_t sgof;
    typedef smart_ptr <scal_3p <strategy_t>, true>	sp_scal_3p_t;
    sp_scal_3p_t scal_3p (params.scal_3p);
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    // Read table for each of region
    for (index_t i = 0; i < idata->sat_region; ++i)
      {
        if ((len = reader->read_table (keyword, sgof, 4)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        scal_3p->get_gas_data()->add_spof (sgof, false);
        sgof.clear ();
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  template <class strategy_t>
  void keyword_manager<strategy_t>::SWFN_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF

    typename strategy_t::item_array_t swfn;
    typedef smart_ptr <scal_3p <strategy_t>, true>	sp_scal_3p_t;
    sp_scal_3p_t scal_3p (params.scal_3p);
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    // Read table for each of region
    for (index_t i = 0; i < idata->sat_region; ++i)
      {
        if ((len = reader->read_table (keyword, swfn, 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        scal_3p->get_water_data()->add_spfn (swfn, i, true);
        swfn.clear ();
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  template <class strategy_t>
  void keyword_manager<strategy_t>::SGFN_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    typename strategy_t::item_array_t sgfn;
    typedef smart_ptr <scal_3p <strategy_t>, true>	sp_scal_3p_t;
    sp_scal_3p_t scal_3p (params.scal_3p);
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    // Read table for each of region
    for (index_t i = 0; i < idata->sat_region; ++i)
      {
        if ((len = reader->read_table (keyword, sgfn, 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        scal_3p->get_gas_data()->add_spfn (sgfn, i, false);
        sgfn.clear ();
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }


  template <class strategy_t>
  void keyword_manager<strategy_t>::SOF3_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    typename strategy_t::item_array_t sof3;
    typedef smart_ptr <scal_3p <strategy_t>, true>	sp_scal_3p_t;
    sp_scal_3p_t scal_3p (params.scal_3p);
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    // Read table for each of region
    for (index_t i = 0; i < idata->sat_region; ++i)
      {
        if ((len = reader->read_table (keyword, sof3, 3)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid argument for keyword %s")
              % reader->get_prefix () % keyword);
          }

        scal_3p->get_water_data()->add_sof3 (sof3, i, true);
        scal_3p->get_gas_data()->add_sof3 (sof3, i, false);
        sof3.clear ();
      }
    BOSOUT (section::read_data, level::medium) <<  keyword << bs_end;
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::SOF2_handler(const std::string &/*keyword*/, keyword_params_t &/*params*/)
  {
    bs_throw_exception ("SOF2 not implemented yet!");
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::EQUIL_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    double dbuf[DOUB_BUF_LEN] = {0};
    char *strt = 0, *end_ptr = 0;
    sp_idata_t idata (params.data, bs_dynamic_cast ());

    //! EQUIL keyword enumerated params

    int i,ii;

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
    if (!idata->rock.size())
      {
        bs_throw_exception (boost::format ("Error in %s: rock properties have not been initialized yet (keyword: %s")
          % reader->get_prefix () % keyword);
      }
    /*if (!(idata->swof.size() && data->sgof.size()) &&
        !(idata->swfn.size() && idata->sgfn.size() && idata->sof3.size()))
    {
      out_s << "Error in " << reader->get_prefix() <<
        ": SCAL tables has not been initialized yet in keyword "
        << keyword;
      KH_ASSERT_EXCEPTION
    }*/
    if (!idata->equil.size())
      {
        bs_throw_exception (boost::format ("Error in %s: EQUIL table has not been initialized yet (keyword: %s")
          % reader->get_prefix () % keyword);
      }

    if (idata->eql_region == 0)
      idata->eql_region = 1;

    for (i = 0; i < (int) idata->eql_region; ++i)
      {
        reader->read_line (buf, CHAR_BUF_LEN);
        if (buf[0] == '/')
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }
        unwrap (buf);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": not valid string format: " << buf;
        //    KH_ASSERT_EXCEPTION
        //  }

        // Values by default
        dbuf[EQUIL_DAT_DEPTH] = -1.0; //Ddat
        dbuf[EQUIL_DAT_PRESS] = 0.0;  //Pdat
        dbuf[EQUIL_WOC_DEPTH] = -1.0; //Dwoc
        dbuf[EQUIL_WOC_PRESS] = 0.0;  //Pwoc
        dbuf[EQUIL_GOC_DEPTH] = 0.0;  // Dgoc
        dbuf[EQUIL_GOC_PRESS] = 0.0;  // Pgoc
        dbuf[EQUIL_RS_TYPE] = 0;      // RS type initialization
        dbuf[EQUIL_RV_TYPE] = 0;      // RV type initialization
        dbuf[EQUIL_NUM_SEC] = 100;      // number of points
        dbuf[EQUIL_COMP_TYPE] = 0;      // composition type initialization
        dbuf[EQUIL_COMP_ARG] = 0;      // composition argument initialization

        strt = buf;
        scanf_d (strt, &end_ptr, &dbuf[EQUIL_DAT_DEPTH]);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": can't read datum depth " << strt;
        //    KH_ASSERT_EXCEPTION
        //  }
        if (dbuf[EQUIL_DAT_DEPTH] < 0) // Ddat
          {
            bs_throw_exception (boost::format ("Error in %s: incorrect value of datum depth %f")
              % reader->get_prefix () % dbuf[EQUIL_DAT_DEPTH]);
          }
        strt = end_ptr;
        scanf_d (strt, &end_ptr, &dbuf[EQUIL_DAT_PRESS]);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": can't read pressure at datum depth " << strt;
        //    KH_ASSERT_EXCEPTION
        //  }
        if (dbuf[EQUIL_DAT_PRESS] < 0) //Pdat
          {
            bs_throw_exception (boost::format ("Error in %s: incorrect pressure value at datum depth %f")
              % reader->get_prefix () % dbuf[EQUIL_DAT_PRESS]);
          }
        strt = end_ptr;
        scanf_d (strt, &end_ptr, &dbuf[EQUIL_WOC_DEPTH]);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": can't read WOC depth " << strt;
        //    KH_ASSERT_EXCEPTION
        //  }
        if (dbuf[EQUIL_WOC_DEPTH] < 0)
          {
            bs_throw_exception (boost::format ("Error in %s: incorrect value of WOC depth %f")
              % reader->get_prefix () % dbuf[EQUIL_WOC_DEPTH]);
          }
        strt = end_ptr;
        scanf_d (strt, &end_ptr, &dbuf[EQUIL_WOC_PRESS]);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": can't read pressure at WOC depth " << strt;
        //    KH_ASSERT_EXCEPTION
        //  }
        strt = end_ptr;
        scanf_d (strt, &end_ptr, &dbuf[EQUIL_GOC_DEPTH]);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": can't read GOC depth " << strt;
        //    KH_ASSERT_EXCEPTION
        //  }
        if (dbuf[EQUIL_GOC_DEPTH] < 0)
          {
            //                rep->print (LOG_READ_SECTION, LOG_ERR,
            //                           GET_TEXT ("Error in %s: incorrect value of GOC depth %lf\n"), r.get_prefix (), dbuf[4]);
            //                  return -8;
            dbuf[EQUIL_GOC_DEPTH] = 0;
          }
        strt = end_ptr;
        scanf_d (strt, &end_ptr, &dbuf[EQUIL_GOC_PRESS]);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": can't read pressure at GOC depth " << strt;
        //    KH_ASSERT_EXCEPTION
        //  }
        strt = end_ptr;
        scanf_d (strt, &end_ptr, &dbuf[EQUIL_RS_TYPE]);
        //if (reader->)
        //  {
        //    out_s << "Error in " << reader->get_prefix() <<
        //    ": can't read RS initialization type at GOC depth " << strt;
        //    KH_ASSERT_EXCEPTION
        //  }
#if 1
        for (ii = 0; ii < EQUIL_TOTAL; ++ii)
          idata->equil[EQUIL_TOTAL * i + ii] = dbuf[ii];
#else
        idata->equil[EQUIL_TOTAL * i + EQUIL_DAT_DEPTH] = dbuf[0]; //Ddat
        idata->equil[EQUIL_TOTAL * i + EQUIL_DAT_PRESS] = dbuf[1]; //Pdat
        idata->equil[EQUIL_TOTAL * i + EQUIL_WOC_DEPTH] = dbuf[2]; //Dwoc
        idata->equil[EQUIL_TOTAL * i + EQUIL_WOC_PRESS] = dbuf[3]; //Pwoc
        idata->equil[EQUIL_TOTAL * i + EQUIL_GOC_DEPTH] = dbuf[4]; //Dgoc
        idata->equil[EQUIL_TOTAL * i + EQUIL_GOC_PRESS] = dbuf[5]; //Pgoc
#endif
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::PRVD_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (idata->eql_region == 0)
      idata->eql_region = 1;

    // allocate memory for prvd
    idata->prvd.resize(idata->eql_region);

    // Read table for each of region
    shared_vector <double> dbuf;
    for (index_t i = 0; i < idata->eql_region; ++i)
      {
        if ((len =reader->read_table (keyword, dbuf, 2)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }

        idata->prvd[i].set_table_len (index_t (len));

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

  template <class strategy_t>
  void keyword_manager<strategy_t>::RSVD_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (idata->eql_region == 0)
      idata->eql_region = 1;

    // allocate memory for rsvd
    idata->rsvd.resize(idata->eql_region);

    // Read table for each of region
    shared_vector <double> dbuf;
    for (index_t i = 0; i < idata->eql_region; ++i)
      {
        if ((len =reader->read_table (keyword, dbuf, 2)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }

        idata->rsvd[i].set_table_len (index_t (len));

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


  template <class strategy_t>
  void keyword_manager<strategy_t>::PBVD_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    sp_idata_t idata (params.data, bs_dynamic_cast ());
    
    if (idata->eql_region == 0)
      idata->eql_region = 1;

    // allocate memory for rsvd
    idata->pbvd.resize(idata->eql_region);

    // Read table for each of region
    shared_vector <double> dbuf;
    for (index_t i = 0; i < idata->eql_region; ++i)
      {
        if ((len =reader->read_table (keyword, dbuf, 2)) < 1)
          {
            bs_throw_exception (boost::format ("Error in %s: not enough valid arguments for keyword %s")
              % reader->get_prefix () % keyword);
          }

        idata->pbvd[i].set_table_len (index_t (len));

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


  template <class strategy_t>
  void keyword_manager<strategy_t>::START_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    sp_this_t km (params.km, bs_dynamic_cast ());
    
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


  template <class strategy_t>
  void keyword_manager<strategy_t>::DATE_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    char buf[CHAR_BUF_LEN] = {0};
    typedef event_manager <strategy_t> event_manager_t;
    typedef smart_ptr <event_manager_t, true>	sp_event_manager_t;
    sp_event_manager_t em (params.em, bs_dynamic_cast ());
    sp_this_t km (params.km, bs_dynamic_cast ());
    
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

  template <class strategy_t>
  void keyword_manager<strategy_t>::DATES_handler(const std::string &keyword, keyword_params_t &params)
  {
    DATE_handler(keyword, params);
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::TSTEP_handler(const std::string &keyword, keyword_params_t &params)
  {
    KH_READER_DEF
    shared_vector <double> dtmp;
    dtmp.resize(MAX_TIME_STEPS_DEF);
    typedef event_manager <strategy_t> event_manager_t;
    typedef smart_ptr <event_manager_t, true>	sp_event_manager_t;
    sp_event_manager_t em (params.em, bs_dynamic_cast ());
    sp_this_t km (params.km, bs_dynamic_cast ());

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

  template <class strategy_t>
  void keyword_manager<strategy_t>::TSTEPS_handler(const std::string &keyword, keyword_params_t &params)
  {
    TSTEP_handler(keyword, params);
  }

//bs stuff
#define KH_SPEC(NAME)\
  template void keyword_manager <base_strategy_di>::NAME  (const std::string &keyword, keyword_params_t &params);\
  template void keyword_manager <base_strategy_fi>::NAME  (const std::string &keyword, keyword_params_t &params);\
  template void keyword_manager <base_strategy_mixi>::NAME  (const std::string &keyword, keyword_params_t &params);

  template keyword_manager<base_strategy_di>::~keyword_manager ();
  template keyword_manager<base_strategy_fi>::~keyword_manager ();
  template keyword_manager<base_strategy_mixi>::~keyword_manager ();

  KH_SPEC(int_array_handler)
  KH_SPEC(float_array_handler)
  KH_SPEC(event_handler)
  KH_SPEC(TITLE_handler)
  KH_SPEC(OIL_handler)
  KH_SPEC(WATER_handler)
  KH_SPEC(GAS_handler)
  KH_SPEC(PROCESS_PARAMS_handler)
  KH_SPEC(RESTART_handler)
  KH_SPEC(REPORTS_handler)
  KH_SPEC(REPORTFILE_handler)
  KH_SPEC(REPORTSCREEN_handler)
  KH_SPEC(STONE1_handler)
  KH_SPEC(STONE2_handler)
  KH_SPEC(RELATIVE_PERM_DEFAULT_handler)
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
  KH_SPEC(WELLDIMS_handler)

  //!TODO: kill next string after debug
  //template class keyword_manager <base_strategy_di>;
  //template class keyword_manager <base_strategy_fi>;

}//nsbs
