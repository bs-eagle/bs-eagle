/** 
 * @file gis_iface.h
 * @brief WELL GIS interface
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-26
 */
#ifndef GIS_IFACE_CT2B01R1

#define GIS_IFACE_CT2B01R1

#include <string>

#include "bs_object_base.h"
#include "conf.h"
#include "table_iface.h"
#include "prop_iface.h"

namespace blue_sky
{
class BS_API_PLUGIN gis_iface : public objbase
  {
    public:
      typedef BS_SP (table_iface)                   sp_table_t;
      typedef BS_SP (prop_iface)                    sp_prop_t;
      typedef BS_SP (gis_iface)                     sp_gis_t;

    public:
      /** 
       * @brief destructor
       */
      virtual ~gis_iface ()
        {}
      /** 
       * @brief return SP to the table
       */
      virtual sp_table_t get_table () = 0;

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop () = 0;

      /** 
       * @brief read data from LAS file
       * 
       * @param fname -- <INPUT> path to the LAS file
       * 
       * @return 0 if ok
       */
      virtual int read_from_las_file (const std::string &fname) = 0;

      virtual sp_gis_t check_serial () const = 0;
      /** 
       * @brief pack(serialize) all information of class to text string 
       * 
       * @return string
       */
      virtual std::string to_str () const = 0; 

      /** 
       * @brief Reastore all class information from input text string
       * 
       * @param s -- <INPUT> string
       */
      virtual void from_str (const std::string &s) = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const = 0;

#endif //BSPY_EXPORTING_PLUGIN
};

}  // end of bluesky name space


#endif /* end of include guard: GIS_IFACE_CT2B01R1 */
