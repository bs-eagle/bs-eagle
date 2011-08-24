/**
* @file keywords.h
* @brief keyword registration
* @author Morozov Andrey
* @date 2008-07-02
*/
#include "bs_bos_core_data_storage_stdafx.h"

/*
#include "event_base.h"
#include "event_manager.h"
*/
#include "keyword_manager.h"
#include "keyword_info_base.h"

namespace blue_sky
  {
  //!TODO: may be use one variable instead of handler_blahblahblah
#define REG_KEYWORD_EXT(KWRD,HNDL) \
    register_keywyfgord(#KWRD, keyword_handler (&this_t::HNDL));

#define REG_KEYWORD(KWRD) \
    register_keyword(#KWRD, keyword_handler (&this_t::KWRD##_handler, 0));

#define REG_ARRAY_KEYWORD(KWRD,T) \
    register_keyword(#KWRD, keyword_handler (&this_t::T##_array_handler, KWRD));

#define REG_INT_ARRAY_KEYWORD(KWRD) \
    REG_ARRAY_KEYWORD(KWRD,int);

#define REG_FLOAT_ARRAY_KEYWORD(KWRD) \
    REG_ARRAY_KEYWORD(KWRD,float);



  
  void keyword_manager::register_keywords()//const sp_event_manager_t &em)
  {
    //Keywords
    REG_KEYWORD(TITLE);
    REG_KEYWORD(OIL);
    REG_KEYWORD(WATER);
    REG_KEYWORD(GAS);
    //REG_KEYWORD(PROCESS_PARAMS);
    REG_KEYWORD(RESTART);
    REG_KEYWORD(REPORTS);
    REG_KEYWORD(REPORTFILE);
    REG_KEYWORD(REPORTSCREEN);
    //REG_KEYWORD(ARITHMETIC);
    REG_KEYWORD(STONE1);
    REG_KEYWORD(STONE2);
    REG_KEYWORD(RELATIVE_PERM_DEFAULT);
    REG_KEYWORD(SCALECRS);
    REG_KEYWORD(UNITS);
    REG_KEYWORD(ROCKCOMP);
    REG_KEYWORD(REGDIMS);
    REG_KEYWORD(REGNUM);
    REG_KEYWORD(EQLDIMS);
    REG_KEYWORD(TABDIMS);
    //REG_KEYWORD(DENSITY);
    REG_KEYWORD(ROCKTAB);
    REG_KEYWORD(START);
    //REG_KEYWORD(PVTO);
    //REG_KEYWORD(PVDO);
    //REG_KEYWORD(PVTW);
    //REG_KEYWORD(PVDG);
    //REG_KEYWORD(ROCK);
/*
    REG_KEYWORD(EQUIL);
    REG_KEYWORD(PRVD);
    REG_KEYWORD(RSVD);
    REG_KEYWORD(PBVD);
    REG_KEYWORD(START);
    REG_KEYWORD(DATE);
    REG_KEYWORD(DATES);
    REG_KEYWORD(TSTEP);
    REG_KEYWORD(TSTEPS);
*/    
    REG_KEYWORD(WELLDIMS);


    // FIXME: npy_intp
    // FIXME: move registration to data_class
    npy_intp dimens[] = {1, 0, 1, 0, 1, 0};
    register_i_pool_keyword ("MPFANUM", dimens, 1, 0);
    register_i_pool_keyword ("EQLNUM", dimens, 1, 0);
    register_i_pool_keyword ("SATNUM", dimens, 1, 0);
    register_i_pool_keyword ("PVTNUM", dimens, 1, 0);
    register_i_pool_keyword ("FIPNUM", dimens, 1, 0);
    register_i_pool_keyword ("BNDNUM", dimens, 2, 0);
    register_i_pool_keyword ("EOSNUM", dimens, 1, 0);
    register_i_pool_keyword ("ROCKNUM", dimens, 1, 0);

    //Pool array
    /*
    //double
    REG_FLOAT_ARRAY_KEYWORD(SWATINIT);
    */

    //BS_ASSERT (em);

    /*
    const type_descriptor &event_td = event_base ::bs_type ();
    const std::vector <type_tuple> &types = BS_KERNEL.registered_types ();
    for (size_t i = 0, cnt = types.size (); i < cnt; ++i)
      {
        const type_descriptor &td = types[i].td_;
        type_descriptor parent_td = td.parent_td ();
        while (!parent_td.is_nil ())
          {
            if (parent_td.stype_ == event_td.stype_)
              {
                register_keyword (td.short_descr_, keyword_handler (&this_t::event_handler));
                break;
              }

            parent_td = parent_td.parent_td ();
          }
      }
    */
  }
  
  
  void keyword_manager::register_plugin_keywords()
  {
    const std::vector <type_tuple> &types = BS_KERNEL.registered_types ();
    sp_keyword_info_base_t keywords;
    const type_descriptor &keyword_info_td = keyword_info_base ::bs_type ();
    sp_objbase keyword_manager (this);
    
    for (size_t i = 0, cnt = types.size (); i < cnt; ++i)
      {
        const type_descriptor &td = types[i].td_;
        type_descriptor parent_td = td;
        while (!parent_td.is_nil ())
          {
            if (parent_td.stype_ == keyword_info_td.stype_)
              {
                keywords = BS_KERNEL.create_object (td);
                BOSOUT (section::keywords, level::low) << boost::format ("Loading keywords from [%s]...") % td.short_descr_ << bs_end;
                keywords->register_keywords (keyword_manager);

                break;
              }

            parent_td = parent_td.parent_td ();
          }
      }
  }

  //!TODO: kill next string after debug
  /*
  template void keyword_manager <base_strategy_did>::register_keywords();//const sp_event_manager_t &em);
  template void keyword_manager <base_strategy_fif>::register_keywords();//const sp_event_manager_t &em);
  template void keyword_manager <base_strategy_dif>::register_keywords();//const sp_event_manager_t &em);
  template void keyword_manager <base_strategy_dld>::register_keywords();//const sp_event_manager_t &em);
  template void keyword_manager <base_strategy_flf>::register_keywords();//const sp_event_manager_t &em);
  template void keyword_manager <base_strategy_dlf>::register_keywords();//const sp_event_manager_t &em);
  
  template void keyword_manager <base_strategy_did>::register_plugin_keywords();
  template void keyword_manager <base_strategy_fif>::register_plugin_keywords();
  template void keyword_manager <base_strategy_dif>::register_plugin_keywords();
  template void keyword_manager <base_strategy_dld>::register_plugin_keywords();
  template void keyword_manager <base_strategy_flf>::register_plugin_keywords();
  template void keyword_manager <base_strategy_dlf>::register_plugin_keywords();
  */
}//ns_bs
