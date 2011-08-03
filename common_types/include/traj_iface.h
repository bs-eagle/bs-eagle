/** 
 * @file traj_iface.h
 * @brief 
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-28
 */
#ifndef TRAJ_IFACE_TCW2C7GV

#define TRAJ_IFACE_TCW2C7GV


#include <string>
#include <sstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "bs_object_base.h"
#include "conf.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()

#include "table_iface.h"

namespace blue_sky
{
class traj_iface : public objbase
  {
    public:
      typedef BS_SP (table_iface)                     sp_table_t;
      typedef boost::archive::text_iarchive           tia_t;
      typedef boost::archive::text_oarchive           toa_t;

    public:
      /** 
       * @brief destructor
       */
      virtual ~traj_iface ()
        {}
      /** 
       * @brief return SP to the table
       */
      virtual sp_table_t get_table () = 0;

      /** 
       * @brief read data from LAS file
       * 
       * @param fname -- <INPUT> path to the LAS file
       * 
       * @return 0 if ok
       */
      virtual int read_from_dev_file (const std::string &fname) = 0;

      virtual void save (toa_t &ar) const = 0;
      virtual void load (tia_t &ar) = 0;
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const = 0;

#endif //BSPY_EXPORTING_PLUGIN
  };
}
#endif /* end of include guard: TRAJ_IFACE_TCW2C7GV */
