/** 
 * @file gis.h
 * @brief WELL GIS storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-26
 */
#ifndef GIS_2XEDKE05

#define GIS_2XEDKE05


#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "gis_iface.h"

namespace blue_sky
{
  
  class BS_API_PLUGIN gis : public gis_iface
    {
    
    public: 

      typedef BS_SP (table_iface)                     sp_table_iface;
      typedef BS_SP (prop_iface)                      sp_prop_iface;
      typedef boost::archive::text_iarchive           tia_t;
      typedef boost::archive::text_oarchive           toa_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~gis ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
      /** 
       * @brief return SP to the table
       */
      virtual sp_table_iface get_table ()
        {
          return sp_table;
        }

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_iface get_prop ()
        {
          return sp_prop;
        }

      /** 
       * @brief read data from LAS file
       * 
       * @param fname -- <INPUT> path to the LAS file
       * 
       * @return 0 if ok
       */
      virtual int read_from_las_file (const std::string &fname);

      virtual void save (toa_t &ar) const
        {
          sp_table->save (ar);
          sp_prop->save (ar);
        }
      virtual void load (tia_t &ar)
        {
          sp_table->load (ar);
          sp_prop->load (ar);
        }
      virtual sp_gis_t check_serial () const;
    public:
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief pack(serialize) all information of class to text string 
       * 
       * @return string
       */
      virtual std::string to_str () const; 

      /** 
       * @brief Reastore all class information from input text string
       * 
       * @param s -- <INPUT> string
       */
      virtual void from_str (const std::string &s);
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;

#endif //BSPY_EXPORTING_PLUGIN
      
    protected:

      int read_ver_info (sp_prop_iface prop, std::string &s);
      int read_wel_info (sp_prop_iface prop, std::string &s, double ver);
      int read_par_info (sp_prop_iface prop, std::string &s, double ver);
      int read_cur_info (sp_prop_iface prop, std::string &s, int n);
      int read_asc_info (std::vector<t_double> &v, std::string &s, int n, std::fstream &file);
      //int read_ver_info (sp_prop_iface prop, const string &s);
      //int read_ver_info (sp_prop_iface prop, const string &s);
      //int read_ver_info (sp_prop_iface prop, const string &s);
      //int read_ver_info (sp_prop_iface prop, const string &s);
      //int read_ver_info (sp_prop_iface prop, const string &s);
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      sp_table_iface sp_table;      //!< table pointer
      sp_prop_iface sp_prop;        //!< ptoperties pointer

      BLUE_SKY_TYPE_DECL (gis);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: GIS_2XEDKE05 */
