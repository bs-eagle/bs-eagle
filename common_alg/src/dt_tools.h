/** 
 * @file dt_tools.h
 * @brief date time tools
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2012-03-01
 */
#ifndef DT_TOOLS_ZA9NDXPL

#define DT_TOOLS_ZA9NDXPL

#include "dt_tools_iface.h"

namespace blue_sky
{

  class BS_API_PLUGIN dt_tools : public dt_tools_iface
    {
    // ------------------------------------
    // METHODS
    // ------------------------------------
    public:
      // destructor
      virtual ~dt_tools ();
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------

      /** 
       * @brief convert C string date to double date
       * 
       * @param buf     -- <INPUT> C string
       * @param date    -- <OUTPUT> date
       * 
       * @return 0 if success
       */
      virtual int cstr2d (const char *buf, double &date) const;

      /** 
       * @brief convert double date to C string in dd.mm.yyyy format
       * 
       * @param date    -- <INPUT> double date
       * @param buf     -- <OUTPUT> C string
       */
      virtual void d2cstr (const double date, char *buf) const;

      /** 
       * @brief convert C string time to double time
       * 
       * @param buf     -- <INPUT> C string
       * @param tm      -- <OUTPUT> double time
       * 
       * @return 0 if success
       */
      virtual int cstr2t (const char *buf, double &tm) const;

      /** 
       * @brief convert double time to C string in hh::mm:ss.ssss format
       * 
       * @param tm      -- <INPUT> double time
       * @param buf     -- <OUTPUT> C string
       */
      virtual void t2cstr (const double tm, char *buf) const;

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
      virtual double d2ymd (const double date, int &year, int &month, int &day) const;

      /** 
       * @brief pack year, month, day to double date
       * 
       * @param year    -- <INPUT> year
       * @param month   -- <INPUT> month
       * @param day     -- <INPUT> day
       * 
       * @return dpouble date
       */
      virtual double ymd2d (const int year, const int month, const int day) const;

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
      virtual double t2hms (const double tm, int &h, int &m, int &s) const;

      /** 
       * @brief pack hour, minute, seconds to double time
       * 
       * @param h       -- <INPUT> hour
       * @param m       -- <INPUT> minute
       * @param s       -- <INPUT> seconds
       * 
       * @return double time
       */
      virtual double hms2t (const int h, const int m, const int s) const;

      /** 
       * @brief convert double date to ECLIPSE C string format
       * 
       * @param date        -- <INPUT> double date
       * @param buf         -- <OUTPUT> C String
       */
      virtual void d2ecl (const double date, char *buf) const;

      /** 
       * @brief convert ECLIPSE C string to double date
       * 
       * @param buf     -- <INPUT> C string
       * 
       * @return double date, if < 0 error occur
       */
      virtual double ecl2d (const char *buf) const;


    public:
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief convert double date to python list [year, month, day, hour, minute, second]
       * 
       * @param d       -- <INPUT> double date
       * 
       * @return python list
       */
      virtual boost::python::list d2date (const double d) const;

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
      virtual double date2d (int year, int month, int day, int hour, int minute, int second) const;

      /** 
       * @brief convert double date to string in dd.mm.yyyy format
       * 
       * @param date    -- <INPUT> double date
       *
       * @return string
       */
      virtual std::string d2str (double d) const;

      /** 
       * @brief convert double time to string in hh:mm:ss.ssss format
       * 
       * @param t       -- <INPUT> double time
       *
       * @return string
       */
      virtual std::string t2str (double t) const;

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

      BLUE_SKY_TYPE_DECL (dt_tools);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: DT_TOOLS_ZA9NDXPL */

