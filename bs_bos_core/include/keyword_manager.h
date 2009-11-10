 /**
 * @file keywords.h
 * @brief declaration of keyword hadling class
 * @author Morozov Andrey
 * @date 2008-07-02
 * */
#ifndef KEYWORD_MANAGER_H_
#define KEYWORD_MANAGER_H_

#include BS_FORCE_PLUGIN_IMPORT ()
//#include "scal_3p.h"
#include "main_def.h"
#include "read_class.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "keyword_manager_iface.h"
#include "date_sim.h"

namespace blue_sky
{

  template <class strategy_t>
  class BS_API_PLUGIN keyword_info_base;

  //! keyword_manager - class-factory which contain a set of handlers for different keywords
  template <class strategy_t>
  class BS_API_PLUGIN keyword_manager: public keyword_manager_iface <strategy_t>
    {
    public:
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
      typedef typename strategy_t::index_t            index_t;
      typedef typename strategy_t::item_t             item_t;
      typedef typename strategy_t::index_array_t      index_array_t;
      typedef typename strategy_t::item_array_t       item_array_t;
      typedef strategy_t                              strategy_type;

      typedef keyword_manager <strategy_t>		        this_t;               //<! self type
      typedef keyword_manager_iface <strategy_t>      base_t;

      typedef idata             							        idata_t;
      typedef rs_mesh_iface <strategy_t>              mesh_iface_t;
      typedef rs_smesh_iface <strategy_t>             smesh_iface_t;
      typedef keyword_info_base <strategy_t>          keyword_info_base_t;
      typedef keyword_params <strategy_t>             keyword_params_t;

      typedef smart_ptr <this_t, true>							  sp_this_t;            //<! smart pointer to self
      typedef smart_ptr <FRead, true>							    sp_reader_t;            //<! smart pointer to reader
      typedef smart_ptr <idata_t, true>						    sp_idata_t;             //<! smart pointr to initial data storage
      typedef smart_ptr <mesh_iface_t, true >         sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true >        sp_smesh_iface_t;
      typedef smart_ptr <keyword_info_base_t, true>   sp_keyword_info_base_t;
      typedef smart_ptr <objbase, true>               sp_objbase;


      typedef typename base_t::handler_t              handler_t;
      typedef typename base_t::keyword_handler        keyword_handler;
      typedef typename base_t::shared_handler_t       shared_handler_t;



      //! typedef for map of handlers, keyword_params &params
      typedef std::map <std::string, keyword_handler> handlers_t;
      typedef std::map <std::string, std::list<std::string> > supported_keywords_t;


      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
      handlers_t handlers; //<!handlers map for active keywords
      supported_keywords_t supported_keywords; //<! all keywords supported by different plugins

      //! temp variables for date and event_handlers
      boost::posix_time::ptime current_date;
      boost::posix_time::ptime starting_date;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      BLUE_SKY_TYPE_DECL_T(keyword_manager);
      ~keyword_manager ();


    public:
      //! initialization used in py constructor
      void init()//const sp_event_manager_t & em)
      {
        this->register_keywords();//em);
        this->register_plugin_keywords();
      }

      //! handle keyword, called from data_manager.read_keyword_file
      void handle_keyword (const std::string &keyword, keyword_params_t &params);

      //! registration of active keyword in factory
      void register_keyword(const std::string &keyword, keyword_handler handler);

      //! registration of active keyword in factory
      void register_keyword (const std::string &keyword, const shared_handler_t &handler);

      //! registration of supported keywords in factory
      void register_supported_keyword(const std::string &keyword, const std::string &provider);

      //! registration of known keywords
      void register_keywords(); //const sp_event_manager_t &em);

      //! registration of known keywords
      void register_plugin_keywords();

      //! return starting date
      boost::posix_time::ptime get_starting_date () {return starting_date;}

      bool
      is_keyword_supported (const std::string &keyword, keyword_params_t &params) const;

      bool 
      is_keyword_activated (const std::string &keyword, keyword_params_t &params) const;

    private:
      //! Following functions are keyword handlers
      //! General functions for pooled keyword handle
      static void int_array_handler                (const std::string &keyword, keyword_params_t &params);
      static void float_array_handler              (const std::string &keyword, keyword_params_t &params);
      //! Handling of event keywords
      static void event_handler                    (const std::string &keyword, keyword_params_t &params);
      //! Named keywords
      static void TITLE_handler                    (const std::string &keyword, keyword_params_t &params);
      static void OIL_handler                      (const std::string &keyword, keyword_params_t &params);
      static void WATER_handler                    (const std::string &keyword, keyword_params_t &params);
      static void GAS_handler                      (const std::string &keyword, keyword_params_t &params);
      static void PROCESS_PARAMS_handler           (const std::string &keyword, keyword_params_t &params);
      static void RESTART_handler                  (const std::string &/*keyword*/, keyword_params_t &/*params*/)
      {
        BOSOUT << "RESTART: NOT_IMPL_YET" << bs_end;
      }
      static void REPORTS_handler                  (const std::string &keyword, keyword_params_t &params);
      static void REPORTFILE_handler               (const std::string &keyword, keyword_params_t &params);
      static void REPORTSCREEN_handler             (const std::string &keyword, keyword_params_t &params);
      static void ARITHMETIC_handler               (const std::string &keyword, keyword_params_t &params);
      static void STONE1_handler                   (const std::string &keyword, keyword_params_t &params);
      static void STONE2_handler                   (const std::string &keyword, keyword_params_t &params);
      static void RELATIVE_PERM_DEFAULT_handler    (const std::string &keyword, keyword_params_t &params);
      static void UNITS_handler                    (const std::string &keyword, keyword_params_t &params);
      static void DIMENS_handler                   (const std::string &keyword, keyword_params_t &params);
      static void ROCKCOMP_handler                 (const std::string &keyword, keyword_params_t &params);
      static void REGDIMS_handler                  (const std::string &keyword, keyword_params_t &params);
      static void REGNUM_handler                   (const std::string &keyword, keyword_params_t &params);
      static void EQLDIMS_handler                  (const std::string &keyword, keyword_params_t &params);
      static void TABDIMS_handler                  (const std::string &keyword, keyword_params_t &params);
      static void COORD_handler                    (const std::string &keyword, keyword_params_t &params);
      static void ZCORN_handler                    (const std::string &keyword, keyword_params_t &params);
      static void MINPV_handler                    (const std::string &keyword, keyword_params_t &params);
      static void MINSV_handler                    (const std::string &keyword, keyword_params_t &params);
      static void DENSITY_handler                  (const std::string &keyword, keyword_params_t &params);
      static void ROCKTAB_handler                  (const std::string &keyword, keyword_params_t &params);
      static void PVTO_handler                     (const std::string &keyword, keyword_params_t &params);
      static void PVDO_handler                     (const std::string &keyword, keyword_params_t &params);
      static void PVTW_handler                     (const std::string &keyword, keyword_params_t &params);
      static void PVDG_handler                     (const std::string &keyword, keyword_params_t &params);
      static void ROCK_handler                     (const std::string &keyword, keyword_params_t &params);
      static void SWOF_handler                     (const std::string &keyword, keyword_params_t &params);
      static void SGOF_handler                     (const std::string &keyword, keyword_params_t &params);
      static void SWFN_handler                     (const std::string &keyword, keyword_params_t &params);
      static void SGFN_handler                     (const std::string &keyword, keyword_params_t &params);
      static void SOF3_handler                     (const std::string &keyword, keyword_params_t &params);
      static void SOF2_handler                     (const std::string &keyword, keyword_params_t &params);
      static void EQUIL_handler                    (const std::string &keyword, keyword_params_t &params);
      static void PRVD_handler                     (const std::string &keyword, keyword_params_t &params);
      static void RSVD_handler                     (const std::string &keyword, keyword_params_t &params);
      static void PBVD_handler                     (const std::string &keyword, keyword_params_t &params);
      static void START_handler                    (const std::string &keyword, keyword_params_t &params);
      static void DATE_handler                     (const std::string &keyword, keyword_params_t &params);
      static void DATES_handler                    (const std::string &keyword, keyword_params_t &params);
      static void TSTEP_handler                    (const std::string &keyword, keyword_params_t &params);
      static void TSTEPS_handler                   (const std::string &keyword, keyword_params_t &params);
      static void WELLDIMS_handler                 (const std::string &/*keyword*/, keyword_params_t &/*params*/)
      {
        BOSOUT << "WELLDIMS: NOT_IMPL_YET" << bs_end;
      }
    };

}//ns bs

#endif //KEYWORD_MANAGER_H_

