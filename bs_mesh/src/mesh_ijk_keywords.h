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
  template <class strategy_t>
  class BS_API_PLUGIN mesh_ijk_keywords: public smesh_keywords<strategy_t>
    {
      public:
        typedef mesh_ijk_keywords<strategy_t>          this_t;
        typedef smesh_keywords<strategy_t>             base_t;
        typedef typename base_t::index_t               index_t;        
        typedef typename base_t::item_t                item_t;
        typedef typename base_t::sp_objbase            sp_objbase;
        typedef typename base_t::keyword_handler       keyword_handler;
        typedef typename base_t::handler_t             handler_t;
        typedef typename base_t::idata_t               idata_t;
        typedef typename base_t::sp_idata_t            sp_idata_t;
        typedef typename base_t::sp_km_iface_t         sp_km_iface_t;
        typedef typename base_t::keyword_params_t      keyword_params_t;
        
        typedef bs_mesh_ijk <strategy_t>               bs_mesh_ijk_t;
        typedef smart_ptr <bs_mesh_ijk_t, true>				 sp_bs_mesh_ijk_t;
        
        typedef rs_mesh_iface <strategy_t>             rs_mesh_iface_t;
        typedef smart_ptr <rs_mesh_iface_t, true>			 sp_mesh_iface_t;
      
      public:
        //! blue-sky class declaration
        BLUE_SKY_TYPE_DECL(mesh_ijk_keywords);
        
        //! default destructor
        virtual ~mesh_ijk_keywords () {};
        
        //! register active and supported keywords
        void register_keywords (sp_objbase &km, std::string provider) const;
        
        //! activate supported keywords
        static void activate_keywords (sp_objbase &km);
        
        //! main handler instatiates class 
        static void mesh_ijk_handler (const std::string &keyword, keyword_params_t &params);
        
    }; 
  };//namespace blue_sky
#endif // MESH_IJK_KEYS_H
