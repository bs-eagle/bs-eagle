/**
 *       \file  keyword_manager.h
 *      \brief  Declaration of keyword handling class
 *     \author  Morozov Andrey
 *       \date  02.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef KEYWORD_MANAGER_H_
#define KEYWORD_MANAGER_H_

#include BS_FORCE_PLUGIN_IMPORT ()
//#include "scal_3p.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "main_def.h"
#include "read_class.h"
#include "data_class.h"

#include "keyword_manager_iface.h"
#include "date_sim.h"

namespace blue_sky
{

  class BS_API_PLUGIN keyword_info_base;

  /**
   * \class keyword_manager
   * \brief Class-factory which contains a set of handlers for different keywords
   * */
  class BS_API_PLUGIN keyword_manager: public keyword_manager_iface 
    {
    public:
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
      typedef keyword_manager 		        this_t;                 //<! self type
      typedef keyword_manager_iface       base_t;

      typedef idata  							        idata_t;
      /*
      typedef rs_mesh_iface               mesh_iface_t;
      typedef rs_smesh_iface              smesh_iface_t;
      */
      typedef keyword_info_base           keyword_info_base_t;
      typedef keyword_params              keyword_params_t;

      typedef smart_ptr <this_t, true>							  sp_this_t;              //<! smart pointer to self
      typedef smart_ptr <FRead, true>							    sp_reader_t;            //<! smart pointer to reader
      typedef smart_ptr <idata_t, true>						    sp_idata_t;             //<! smart pointr to initial data storage
      /*
      typedef smart_ptr <mesh_iface_t, true >         sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true >        sp_smesh_iface_t;
      */
      typedef smart_ptr <keyword_info_base_t, true>   sp_keyword_info_base_t;
      typedef smart_ptr <objbase, true>               sp_objbase;
      

      typedef base_t::handler_t              handler_t;
      typedef base_t::keyword_handler        keyword_handler;
      typedef base_t::shared_handler_t       shared_handler_t;



      //! typedef for map of handlers, keyword_params &params
      typedef std::map <std::string, keyword_handler> handlers_t;
      typedef std::map <std::string, std::list<std::string> > supported_keywords_t;


    public:
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
      handlers_t handlers;                                //!< Handlers map for active keywords
      supported_keywords_t supported_keywords;            //!< All keywords supported by different plugins

      boost::posix_time::ptime current_date;              //!< temp variable for date and event_handlers
      boost::posix_time::ptime starting_date;             //!< temp variable for date and event_handlers

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL_T(keyword_manager);

      /**
       * \brief  dtor
       * */
      ~keyword_manager ();


    public:
      /**
       * \brief  Registers built-in keyword and keywords from plugins
       * */
      void 
      init()
      {
        //this->register_keywords();
        this->register_plugin_keywords();
      }

      /**
       * \brief  Handles keyword, call proper handler
       * \return 
       * */
      void 
      handle_keyword (const std::string &keyword, keyword_params_t &params);

      /**
       * \brief  Registers active keyword in factory
       * \param  keyword Keyword name
       * \param  handler Keyword handler description
       * */
      void 
      register_keyword (const std::string &keyword, keyword_handler handler);

      /**
       * \brief  Registers active keyword in factory
       * \param  keyword Keyword name
       * \param  handler Instance of object which implements 
       *                 keyword handler interface
       * */
      void 
      register_keyword (const std::string &keyword, const shared_handler_t &handler, bool replace_existing);
      
      //! registration of active integer pool keyword in factory
      void register_i_pool_keyword(const std::string &keyword, int *dimens, t_int def_value, handler_t external_handler = 0);
      
      //! registration of active floating point pool keyword in factory
      void register_fp_pool_keyword(const std::string &keyword, int *dimens, t_float def_value, handler_t external_handler = 0);
      
      //! python registration of active floating point pool keyword in factory
      void py_register_fp_pool_keyword (const std::string keyword, boost::python::list dimens, t_float def_value);

      /**
       * \brief  Registers supported keyword by plugins in factory
       * \param  keyword Keyword name
       * \param  provider
       * */
      void 
      register_supported_keyword(const std::string &keyword, const std::string &provider);

      /**
       * \brief  Registers built-in keywords
       * */
      void 
      register_keywords(); 

      /**
       * \brief  Registers keyword from plugins
       * */
      void 
      register_plugin_keywords();
      
      boost::python::list py_list_active_keywords();
      
      boost::python::list py_list_supported_keywords();

      /**
       * \brief  Returns starting date
       * \return Starting date
       * */
      boost::posix_time::ptime 
      get_starting_date () {return starting_date;}

      /**
       * \brief  Returns true if keyword supported by plugins
       * \param  keyword Keyword name
       * \param  params
       * \return True if keyword supported by plugins
       * */
      bool
      is_keyword_supported (const std::string &keyword, keyword_params_t &params) const;

      /**
       * \brief  Returns true if any handler binded to keyword
       * \param  keyword Keyword name
       * \param  params
       * \return True if any handler binded to keyword
       * */
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
      //static void PROCESS_PARAMS_handler           (const std::string &keyword, keyword_params_t &params);
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
      /*
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
      */
      static void WELLDIMS_handler                 (const std::string &/*keyword*/, keyword_params_t &/*params*/)
      {
        BOSOUT << "WELLDIMS: NOT_IMPL_YET" << bs_end;
      }
    };

}//ns bs

#endif //KEYWORD_MANAGER_H_

