/**
 *       \file  event_base.h
 *      \brief  Base class for model events
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef EVENT_BASE_H_
#define EVENT_BASE_H_
#include "data_class.h"
// WTF??
#include "data_storage_interface.h"

#include BS_FORCE_PLUGIN_IMPORT()
#include "named_pbase_access.h"
#include "rs_mesh_iface.h"
#include BS_STOP_PLUGIN_IMPORT()

namespace blue_sky
  {

  class reservoir;
  class calc_model;
  class idata;

  /**
   * \class event_base
   * \brief Base class for model events
   * */
  class BS_API_PLUGIN event_base : public objbase
    {
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
    public:
      //TODO:change objbase to reservoir
      typedef reservoir          reservoir_t;
      typedef rs_mesh_iface      mesh_iface_t;
      typedef calc_model         calc_model_t;

      typedef smart_ptr <mesh_iface_t , true> sp_mesh_iface_t;
      typedef smart_ptr <reservoir_t, true>   sp_top;
      typedef smart_ptr <calc_model_t, true>  sp_calc_model_t;
      typedef smart_ptr <named_pbase, true>   sp_named_pbase;

      typedef event_base         self_t;
      typedef self_t                          this_t;

      typedef named_pbase::ndt_t              ndt_t;
      typedef named_pbase::named_pbase_idxs   named_pbase_idxs;

      //additional
    public:

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL(event_base);

    public:
      //! destructor
      virtual ~event_base ();

      /**
       * \brief  Applies actions from event
       * \param  top Pointer to reservoir instance
       * \param  mesh Pointer to mesh instance
       * \param  calc_model Pointer to calc_model instance
       * \return
       * */
      virtual void
      apply (const sp_top &top, const sp_mesh_iface_t &mesh,
             const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &data) const;

      /**
       * \brief  Inits event
       * \param  event_name Name of event
       * \param  params String with event params
       * */
      void
      init (const std::string &event_name, const std::string &params);

      /**
       * \brief  Adds subevent
       * \param  event_name Name of event
       * \param  params String with event params
       * */
      void
      add_next_line (const std::string &event_name, const std::string &params);

      /**
       * \brief  Checks is keyword (or event) is multiline
       * \return True is event is multiline
       * */
      virtual bool
      is_multi_line () const
      {
        return false;
      }

    protected:



      /**
       * \brief  Returns pointer to main_params instance
       * \return Instance of named_pbase
       * */
      virtual sp_named_pbase
      main_params () const
      {
        return 0;
      }

      /**
       * \brief  Creates params for new subevent and returns it
       * \return Instance of named_pbase
       * */
      virtual sp_named_pbase
      add_next_line_params ()
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Parses event parameters
       * \param  params_str String with event parameters that should be parsed
       * \param  params Parameters holder
       * */
      void
      parse (const std::string &params_str, const sp_named_pbase &params);

      /**
       * \brief  Saves default value in parameters holder
       * \param  params Parameters holder
       * \param  c _not_used_
       * */
      void
      save_def (const sp_named_pbase &params, const char c = '_'); //!< skip a default value

      /**
       * \brief  Resets num (internal parser state)
       * \param  params Parameters holder
       * \param  begin Begin of chunk of parsed string
       * \param  end End of chunk of parsed string
       * */
      void
      clear_num (const sp_named_pbase &params, const char *begin = 0, const char *end = 0);

      //!<simple ..
      /**
       * \brief  Saves parsed value to holder
       * \param  params Parameters holder
       * \param  first Begin of chunk of parsed string
       * \param  last End of chunk of parsed string
       * */
      void
      save_value (const sp_named_pbase &params, char const* first, char const* last); //< function saves any income value

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    private:
#ifdef _DEBUG
      /**
       * \class debug_data
       * \brief Stores debug data for parser
       * */
      struct debug_data
      {
        std::string               params;                     //!< event parameters
        std::list <std::string>   next_line_params;           //!< next lines params
      };

      debug_data                  debug_data_;                //!< debug data
#endif

      int                         num;                        //!< temporary variables to store a number of same parameters
      int                         index;                      //!< index of parameter which is parsed in this moment
      bool                        no_error_;
    };

  /**
   * \brief  Registers event_base types in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return True if all types registered successfully
   * */
  bool
  event_base_register_types (const plugin_descriptor &pd);

}//namespace blue_sky

#endif//EVENT_BASE_H_

