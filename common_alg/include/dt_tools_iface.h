/**
 * @file dt_tools_iface.h
 * @brief interface for date and time tools
 * @author Oleg Borschuk
 * @version
 * @date 2012-03-01
 */
#ifndef DT_TOOLS_IFACE_28OJZ2EW

#define DT_TOOLS_IFACE_28OJZ2EW




#include <string>

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python/list.hpp>
#endif //BSPY_EXPORTING_PLUGIN

#include "bs_object_base.h"
#include "conf.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
{
class dt_tools_iface : public objbase
  {
    //////////////////////////////////////////
    // METHODS
    // ---------------------------------------

    public:
      /**
       * @brief destructor
       */
      virtual ~dt_tools_iface ()
        {}

      /** 
       * @brief convert C string date to double date
       * 
       * @param buf     -- <INPUT> C string
       * @param date    -- <OUTPUT> date
       * 
       * @return 0 if success
       */
      virtual int cstr2d (const char *buf, double &date) const = 0;

      /** 
       * @brief convert double date to C string in dd.mm.yyyy format
       * 
       * @param date    -- <INPUT> double date
       * @param buf     -- <OUTPUT> C string
       */
      virtual void d2cstr (const double date, char *buf) const = 0;

      /** 
       * @brief convert C string time to double time
       * 
       * @param buf     -- <INPUT> C string
       * @param tm      -- <OUTPUT> double time
       * 
       * @return 0 if success
       */
      virtual int cstr2t (const char *buf, double &tm) const = 0;

      /** 
       * @brief convert double time to C string in hh::mm:ss.ssss format
       * 
       * @param tm      -- <INPUT> double time
       * @param buf     -- <OUTPUT> C string
       */
      virtual void t2cstr (const double tm, char *buf) const = 0;

      /** 
       * @brief unpack double date to year, month, day
       * 
       * @param date    -- <INPUT> double date
       * @param year    -- <OUTPUT> year
       * @param month   -- <OUTPUT> month
       * @param day     -- <OUTPUT> day
       * 
       * @return double time (double date - date)
       */
      virtual double d2ymd (const double date, int &year, int &month, int &day) const = 0;

      /** 
       * @brief pack year, month, day to double date
       * 
       * @param year    -- <INPUT> year
       * @param month   -- <INPUT> month
       * @param day     -- <INPUT> day
       * 
       * @return dpouble date
       */
      virtual double ymd2d (const int year, const int month, const int day) const = 0;

      /** 
       * @brief unpack double time to hour, minute, sec,
       * 
       * @param tm          -- <INPUT> double time
       * @param h           -- <OUTPUT> hour
       * @param m           -- <OUTPUT> minute
       * @param s           -- <OUTPUT> sec
       * 
       * @return mileseconds
       */
      virtual double t2hms (const double tm, int &h, int &m, int &s) const = 0;

      /** 
       * @brief pack hour, minute, seconds to double time
       * 
       * @param h       -- <INPUT> hour
       * @param m       -- <INPUT> minute
       * @param s       -- <INPUT> seconds
       * 
       * @return double time
       */
      virtual double hms2t (const int h, const int m, const int s) const = 0;

      /** 
       * @brief convert double date to ECLIPSE C string format
       * 
       * @param date        -- <INPUT> double date
       * @param buf         -- <OUTPUT> C String
       */
      virtual void d2ecl (const double date, char *buf) const = 0;
      
      /** 
       * @brief convert ECLIPSE C string to double date
       * 
       * @param buf     -- <INPUT> C string
       * 
       * @return double date, if < 0 error occur
       */
      virtual double ecl2d (const char *buf) const = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      
      /** 
       * @brief convert double date to python list [year, month, day, hour, minute, second]
       * 
       * @param d       -- <INPUT> double date
       * 
       * @return python list
       */
      virtual boost::python::list d2date (const double d) const = 0;

      /** 
       * @brief convert [year, month, day, hour, minute, second] to double date
       * 
       * @param year    -- <INPUT> year
       * @param month   -- <INPUT> month
       * @param day     -- <INPUT> day
       * @param hour    -- <INPUT> hour
       * @param minute  -- <INPUT> minute
       * @param second  -- <INPUT> second
       * 
       * @return double date
       */
      virtual double date2d (int year, int month, int day, int hour, int minute, int second) const = 0;

      /** 
       * @brief convert double date to string in dd.mm.yyyy format
       * 
       * @param date    -- <INPUT> double date
       *
       * @return string
       */
      virtual std::string d2str (double d) const = 0;

      /** 
       * @brief convert double time to string in hh:mm:ss.ssss format
       * 
       * @param t       -- <INPUT> double time
       *
       * @return string
       */
      virtual std::string t2str (double t) const = 0;

      /**
       * @brief python print wrapper
       *
       * @return return table description
       */
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN
};

}  // end of bluesky name space

#endif /* end of include guard: DT_TOOLS_IFACE_28OJZ2EW */

