/** 
 * @file pvt_iface.h
 * @brief PVT interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-25
 */
#ifndef PVT_IFACE_ZT72GBDA

#define PVT_IFACE_ZT72GBDA


#include "bs_assert.h"
#include "bs_tree.h"
#include "smart_ptr.h"
#include "bs_array.h"
#include "conf.h"
#include "prop_iface.h"
#include "table_iface.h"

#include <string>

namespace blue_sky
{
  /** 
   * @brief PVT interface
   */
  class pvt_iface: public objbase
  {
    public:
      typedef smart_ptr<prop_iface, true>       sp_prop_t;
      typedef smart_ptr<table_iface, true>      sp_table_t;


      /** 
       * @brief destructor
       */
       virtual ~pvt_iface () {}

       /** 
        * @brief check class consistency 
        * 
        * @return 0 if ok, < 0 error found
        */
       virtual int check () const = 0;

       /** 
        * @brief update calculated part of PVT
        * 
        * @return 0 if ok, <0 if error occur
        */
       virtual int update () = 0;

       /** 
        * @brief return surface density of the phase
        */
       virtual t_double get_surface_density () const = 0;

       /** 
        * @brief return input properties for a given region
        */
       virtual sp_prop_t get_prop () = 0;

       /** 
        * @brief check and set new property
        * 
        * @param new_prop -- <INPUT> given property
        */
       virtual void set_prop (sp_prop_t new_prop) = 0;

       /** 
        * @brief return table to the input data
        * 
        */
       virtual sp_table_t get_table () = 0; 

       /** 
        * @brief check and set given table 
        * 
        * @param new_table -- <INPUT> given table
        */
       virtual void set_table (sp_table_t new_table) = 0;


              
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

  };
};

#endif /* end of include guard: PVT_IFACE_ZT72GBDA */

