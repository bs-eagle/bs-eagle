#ifndef MESH_IJK_KEYS_H
#define MESH_IJK_KEYS_H
/*!
	\file mesh_ijk_keywords.h
  \brief This file declares keyword handlers for bs reservoir simulation mesh IJK
  \author Mark Khait
	\date 2009-08-07
 */

#include "rs_smesh_keywords.h"
#include "bs_mesh_ijk.h"

namespace blue_sky
  {

 //! keyword_register_iface - interface for register plugins keywords
  class BS_API_PLUGIN mesh_ijk_keywords: public smesh_keywords
    {
      public:
        typedef mesh_ijk_keywords          this_t;
        typedef smesh_keywords             base_t;
        
        typedef bs_mesh_ijk                bs_mesh_ijk_t;
        typedef smart_ptr <bs_mesh_ijk_t, true>        sp_bs_mesh_ijk_t;
        
        typedef rs_mesh_iface              rs_mesh_iface_t;
        typedef smart_ptr <rs_mesh_iface_t, true>      sp_mesh_iface_t;
      
      public:
        //! blue-sky class declaration
        BLUE_SKY_TYPE_DECL(mesh_ijk_keywords);
        
        //! default destructor
        virtual ~mesh_ijk_keywords () {};
        
        //! register active and supported keywords
        void register_keywords (sp_objbase &km, std::string provider) const;
        
        //! activate supported keywords
        static void activate_keywords (sp_km_iface_t keyword_manager);
        
        //! main handler instatiates class 
        static void mesh_ijk_reactor (const std::string &keyword, keyword_params_t &params);
        
    }; 
  };//namespace blue_sky
#endif // MESH_IJK_KEYS_H
