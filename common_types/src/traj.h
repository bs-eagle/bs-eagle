/** 
 * @file traj.h
 * @brief wellbore trajectory implementation
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-28
 */
#ifndef TRAJ_DAKO7FWL

#define TRAJ_DAKO7FWL


#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "traj_iface.h"
#include "bs_serialize_decl.h"

namespace blue_sky
{

  class BS_API_PLUGIN traj : public traj_iface
    {

    public:

      typedef BS_SP (table_iface)                   sp_table_t;
      typedef BS_SP (prop_iface)                    sp_prop_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~traj ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
      /** 
       * @brief return SP to the table
       */
      virtual sp_table_t get_table ()
        {
          return sp_table;
        }

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop ()
        {
          return sp_prop;
        }

      /** 
       * @brief read data from DEV file
       * 
       * @param fname -- <INPUT> path to the DEV file
       * 
       * @return 0 if ok
       */
      virtual int read_from_dev_file (const std::string &fname);
    public:
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
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;

#endif //BSPY_EXPORTING_PLUGIN

    protected:
      // ------------------------------
      // VARIABLES
      // ------------------------------
      sp_table_t sp_table;      //!< table pointer
      sp_prop_t  sp_prop;    //!< properties pointer

      BLUE_SKY_TYPE_DECL (traj);
      friend class bs_serialize;
    };

}; //end of blue_sky namespace

#endif /* end of include guard: TRAJ_DAKO7FWL */
