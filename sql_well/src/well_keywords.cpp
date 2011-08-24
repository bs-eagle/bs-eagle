/**
 *       \file  well_keywords.cpp
 *      \brief  keywords for wells
 *     \author  Mark Khait
 *       \date  05.08.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "well_keywords.hpp"
#include "keyword_manager_iface.h"
#include "well_pool_iface.h"
#include "prop_iface.h"



namespace blue_sky 
{
  
  
  well_keywords::well_keywords (bs_type_ctor_param)
  {
  }

  well_keywords::well_keywords (const well_keywords &src)
  : bs_refcounter (src), keyword_info_base (src)
  {
  }

  void
  well_keywords::register_keywords (sp_objbase &km, std::string provider) const
  {
    BS_SP (keyword_manager_iface) keyword_manager (km, bs_dynamic_cast ());
    BS_ASSERT (keyword_manager);
    std::vector<std::string> names(1);
  
    names[0] = "csv_schedule";
  
    keyword_manager->register_prop_keyword ("CSV_SCHEDULE", "S", names, &this_t::CSV_SHEDULE_reactor);
    //keyword_manager->register_keyword ("CSV_SHEDULE", keyword_handler(CSV_SHEDULE_reader, CSV_SHEDULE_reactor);
   
  }
  
  void
  well_keywords::CSV_SHEDULE_reactor (std::string const &keyword, keyword_params &params)
  {
    char prop_name[50];
    char prop_name_last[50];
    t_int j = 1;
    
    sprintf (prop_name, "csv_schedule");
    
    while (1)
      {
        try
          {
            params.hdm->get_prop()->get_s(prop_name);
            sprintf (prop_name_last, "%s", prop_name);
            sprintf (prop_name, "%s_%d", "csv_schedule", j);
            j++; 
          }
        catch (...)
          {
            break;
          }
      }
    std::string csv_schedule = params.hdm->get_prop()->get_s (prop_name_last);
    t_double starting_date = params.hdm->get_prop()->get_f ("starting_date");
    params.hdm->get_well_pool()->read_from_ascii_file(csv_schedule, starting_date);
  }

  BLUE_SKY_TYPE_STD_CREATE (well_keywords);
  BLUE_SKY_TYPE_STD_COPY (well_keywords);
  BLUE_SKY_TYPE_IMPL (well_keywords, keyword_info_base, "Well Keywords", "well_keywords", "well_keywords");

}
