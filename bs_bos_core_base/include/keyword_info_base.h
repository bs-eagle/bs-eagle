 /**
 * @file keywords.h
 * @brief plugin keywords info base class
 * @author Mark Khait
 * @date 2009-08-03
 * */
#ifndef KEYWORD_INFO_BASE_H_
#define KEYWORD_INFO_BASE_H_

#include "keyword_manager_iface.h"


namespace blue_sky
{
  
  //! keyword_register_iface - interface for register plugins keywords
  class BS_API_PLUGIN keyword_info_base: public objbase
    {
      public:
      
        typedef keyword_manager_iface                km_t;
        typedef smart_ptr <objbase, true>            sp_objbase;
        typedef km_t::handler_t             handler_t;
        typedef km_t::keyword_handler       keyword_handler;


      public:
      
        //! blue-sky class declaration
        virtual void register_keywords (sp_objbase &/*km*/, std::string /*provider*/ = "") const {};


        BLUE_SKY_TYPE_DECL(keyword_info_base);
    };
}//ns bs

#endif //KEYWORD_INFO_BASE_H_

