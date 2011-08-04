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

namespace blue_sky
{
  
  class BS_API_PLUGIN traj : public traj_iface
    {
    
    public: 

      typedef BS_SP (table_iface)                     sp_table_t;
      typedef boost::archive::binary_iarchive           tia_t;
      typedef boost::archive::binary_oarchive           toa_t;

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
       * @brief read data from DEV file
       * 
       * @param fname -- <INPUT> path to the DEV file
       * 
       * @return 0 if ok
       */
      virtual int read_from_dev_file (const std::string &fname);
      virtual void save (toa_t &ar) const
        {
          sp_table->save (ar);
        }
      virtual void load (tia_t &ar)
        {
          sp_table->load (ar);
        }
    public:
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
    protected:
      sp_table_t sp_table;      //!< table pointer

      BLUE_SKY_TYPE_DECL (traj);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: TRAJ_DAKO7FWL */
