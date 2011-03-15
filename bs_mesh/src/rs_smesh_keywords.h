#ifndef RS_SMESH_KEYS_H
#define RS_SMESH_KEYS_H
/*!
	\file rs_smesh_keywords.h
  \brief This file declares keyword handlers for bs reservoir simulation structured meshes keywords
  \author Mark Khait
	\date 2009-08-07
 */

#include "keyword_info_base.h"
#include "rs_smesh_iface.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "data_class.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  class BS_API_PLUGIN smesh_keywords: public keyword_info_base
    {
      public:
      
        typedef smesh_keywords             this_t;
        typedef keyword_info_base          base_t;

        typedef keyword_params             keyword_params_t;
        
        typedef rs_mesh_iface              rs_mesh_iface_t;
        typedef smart_ptr <rs_mesh_iface_t, true>      sp_mesh_iface_t;
        
        typedef rs_smesh_iface             rs_smesh_iface_t;
        typedef smart_ptr <rs_smesh_iface_t, true>     sp_smesh_iface_t;
        
        typedef keyword_manager_iface      km_iface_t;
        typedef smart_ptr <km_iface_t, true>           sp_km_iface_t;
        
        typedef idata                      idata_t;
        typedef smart_ptr <idata_t, true>	             sp_idata_t;
        
        typedef FRead                                  reader_t;
        typedef smart_ptr <reader_t, true>	           sp_reader_t;
      
      public:
        //! blue-sky class declaration
        BLUE_SKY_TYPE_DECL(smesh_keywords);
        
        //! default destructor
        virtual ~smesh_keywords () {};
        
        //! register active and supported keywords
        void register_keywords (sp_objbase &km, std::string provider) const;
        
        //! activate supported keywords
        static void activate_keywords (sp_objbase &km);
        
        //! mesh dimensions handler 
        static void DIMENS_handler (const std::string &keyword, keyword_params_t &params);
    }; 
};//namespace blue_sky
#endif // RS_SMESH_KEYS_H
