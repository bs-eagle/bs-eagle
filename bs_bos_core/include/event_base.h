/**
 * @file event_base.h
 * @brief
 * @author Morozov Andrey
 * @date 2008-06-07
 * */

#ifndef EVENT_BASE_H_
#define EVENT_BASE_H_

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"
#include "data_storage_interface.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class reservoir;

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  class BS_API_PLUGIN event_base : public objbase //: public named_pbase//, public combase
    {
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
    public:
      //TODO:change objbase to reservoir
      typedef reservoir <strategy_t>          reservoir_t;
      typedef rs_mesh_iface <strategy_t>      mesh_iface_t;
      typedef calc_model <strategy_t>					calc_model_t;

      typedef smart_ptr <mesh_iface_t , true> sp_mesh_iface_t;
      typedef smart_ptr <reservoir_t, true>   sp_top;
      typedef smart_ptr <calc_model_t, true>	sp_calc_model_t;
      typedef smart_ptr <named_pbase, true>   sp_named_pbase;

      typedef event_base <strategy_t>         self_t;
      typedef self_t                          this_t;

      typedef named_pbase::ndt_t              ndt_t;
      typedef named_pbase::named_pbase_idxs   named_pbase_idxs;

      //additional
    public:

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
    public:
      BLUE_SKY_TYPE_DECL(event_base);

    public:
      //! destructor
      virtual ~event_base ();

      virtual void apply(const sp_top &top, const sp_mesh_iface_t &mesh,
                   const sp_calc_model_t &calc_model) const;										//!< virtual apply function

      //!< initialization
      void 
      init (const std::string &event_name, const std::string &params);
      
      void
      add_next_line (const std::string &event_name, const std::string &params);

      //!< check is keyword (or event) is multiline
      virtual bool 
      is_multi_line () const
      {
        return false;
      }

    protected:
      


      virtual sp_named_pbase
      main_params () const
      {
        return 0;
      }

      virtual sp_named_pbase
      add_next_line_params () 
      {
        bs_throw_exception ("PURE CALL");
      }


      //!< parameters string parsing
      void 
      parse (const std::string &params_str, const sp_named_pbase &params);

      //!< functions used by boost spirit parser
      void
      save_def (const sp_named_pbase &params, const char c = '_'); //!< skip a default value

      //!< clear num
      void
      clear_num (const sp_named_pbase &params, const char *begin = 0, const char *end = 0);
      
      //!<simple ..
      void 
      save_value(const sp_named_pbase &params, char const* first, char const* last); //< function saves any income value

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    private:
#ifdef _DEBUG
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

  bool
  event_base_register_types (const plugin_descriptor &pd);

}//namespace blue_sky

#endif//EVENT_BASE_H_

