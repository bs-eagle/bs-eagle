/** 
 * @file dt_tools.cpp
 * @brief Reader for the BOS ascii files
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2012-03-01
 */

#include "dt_tools.h"
#include "bs_kernel.h"
//#include "localization.h"
#include <math.h>
#include "toupper.h"

using namespace boost;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN

/*! \brief Nonzero if YEAR is a leap year (every 4 years,\n
except every 100th isn't, and every 400th is).  */
#define IS_LEAP_YEAR(YEAR) \
  ((YEAR) % 4 == 0 && ((YEAR) % 100 != 0 || (YEAR) % 400 == 0))


namespace blue_sky
{

  dt_tools::dt_tools (bs_type_ctor_param)
    {
    }
  dt_tools::dt_tools (const dt_tools& /*rhs*/)
        : bs_refcounter ()
    {
    }
  dt_tools::~dt_tools ()
    {
    }

  int 
  dt_tools::cstr2d (const char *buf, double &date) const
    {
      int day, month, year;
      
      date = ecl2d (buf);
      if (date < 0)
        {
          if (sscanf (buf, "%d.%d.%d", &day, &month, &year) != 3)
            {
              fprintf (stderr, "Error: invalid date format: %s\n", buf);
              return -1;
            }
          date = ymd2d (year, month, day);
        }
      return 0;
    }

  void 
  dt_tools::d2cstr (const double date, char *buf) const
    {
      int d, m, y;
      d2ymd (date, y, m, d);
      sprintf (buf, "%02d.%02d.%04d", d, m, y);
    }

  int 
  dt_tools::cstr2t (const char *buf, double &tm) const
    {
      int h, m;
      double s;
      
      if (sscanf (buf, "%d:%d:%lf", &h, &m, &s) != 3)
        {
          fprintf (stderr, "Error: invalid time format: %s\n", buf);
          return -1;
        }
     
      // Check month and year
      
      if (h < 0 || h > 23)
        {
          fprintf (stderr, "Error: hours should be 0..23: %s\n", buf);
          return -1;
        }
      
      if (m < 0 || m > 59)
        {
          fprintf (stderr, "Error: minutes should be 0..59: %s\n", buf);
          return -1;
        } 
      if (s < 0 || s >= 60.0)
        {
          fprintf (stderr, "Error: seconds should be 0..59: %s\n", buf);
          return -1;
        } 
      tm = hms2t (h, m, s); //((double)h + ((double)m + (double)s / 60.0) / 60.0) / 24.0;
      return 0;
    }

  void 
  dt_tools::t2cstr (const double tm, char *buf) const
    {
      int h, m, s;
      //double ss;

      //ss = d2hms (tm, h, m, s);
      t2hms (tm, h, m, s);
      //ss += s;

      sprintf (buf, "%02d:%02d:%02d", h, m, s);
    }

  double 
  dt_tools::d2ymd (const double date, int &year, int &month, int &day) const
    {
      static const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
      double next_date, datetime, d;
      datetime = date;
      d = floor(date);
      for (year = 0; d > 0; year++)
        {
          if (IS_LEAP_YEAR (1900 + year))
            next_date = d - 366;
          else
            next_date = d - 365;
          if (next_date <= 0)
            break;
          d = next_date;
        }

      for (month = 0; d > 0; month++)
        {
          next_date = d - days_in_month[month];
          if (IS_LEAP_YEAR (1900 + year) && month == 1)
            next_date--;
          if (next_date <= 0)
            break;
          d = next_date;
        }
      
      day = (int) d;
      month++;
      year += 1900;
        
      return (datetime - floor(datetime));
    }

  double 
  dt_tools::ymd2d (const int y_, const int m_, const int d_) const
    {
      int days_in_current_month;
      static const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
      int i;
      double date;
      int month = m_, year = y_, day = d_;
      
      // Check month and year
      
      if (year < 1900)
        {
          fprintf (stderr, "Error: years earlier 1900 not supported\n");
          return -1;
        }
      
      if (month < 1 || month > 12)
        {
          fprintf (stderr, "Error: invalid (nonexistent) date\n");
          return -1;
        } 
      
      date = 0;
      
      for (i = 1900; i < year; i++)
        {
        if (IS_LEAP_YEAR (i))
          date += 366;
        else
          date += 365;
        }
        
      month--;  
      
      for (i = 0; i < month; i++)
        {
          date += days_in_month[i];
          if (IS_LEAP_YEAR (year) && i == 1) // February of leap year
            date++;
        }
      date += day;
      
      days_in_current_month = days_in_month[month];
      if (IS_LEAP_YEAR (year) && month == 1) // February of leap year
        days_in_current_month++;
      
      
      // check day
      if (day < 1 || day > days_in_current_month)
        {
          fprintf (stderr, "Error: invalid (nonexistent) date\n");
          return -1;
        }
      
      return date;
    }

  double 
  dt_tools::t2hms (const double tm, int &h, int &m, int &s) const
    {
      double d, n, hh, mm, ss;
      
      d = modf (tm + 1.0e-13, &n);
      d = modf (d * 24.0, &hh);
      d = modf (d * 60.0, &mm);
      d = modf (d * 60.0, &ss);
      h = (int)hh;
      m = (int)mm;
      s = (int)ss;

      return d;
    }

  double 
  dt_tools::hms2t (const int h, const int m, const int s) const
    {
      return ((double)h + ((double)m + (double)s / 60.0) / 60.0) / 24.0;
    }

  void 
  dt_tools::d2ecl (const double date, char *buf) const
    {
      int d, m, y;
      d2ymd (date, y, m, d);
      
      switch (m)
      {
        case 1:
          sprintf (buf, "%2.2d \'JAN\' %4.4d /", d, y);
          break;
        case 2:
          sprintf (buf, "%2.2d \'FEB\' %4.4d /", d, y);
          break;
        case 3:
          sprintf (buf, "%2.2d \'MAR\' %4.4d /", d, y);
          break;
        case 4:
          sprintf (buf, "%2.2d \'APR\' %4.4d /", d, y);
          break;
        case 5:
          sprintf (buf, "%2.2d \'MAY\' %4.4d /", d, y);
          break;
        case 6:
          sprintf (buf, "%2.2d \'JUN\' %4.4d /", d, y);
          break;
        case 7:
          sprintf (buf, "%2.2d \'JUL\' %4.4d /", d, y);
          break;
        case 8:
          sprintf (buf, "%2.2d \'AUG\' %4.4d /", d, y);
          break;
        case 9:
          sprintf (buf, "%2.2d \'SEP\' %4.4d /", d, y);
          break;
        case 10:
          sprintf (buf, "%2.2d \'OCT\' %4.4d /", d, y);
          break;
        case 11:
          sprintf (buf, "%2.2d \'NOV\' %4.4d /", d, y);
          break;
        case 12:
          sprintf (buf, "%2.2d \'DEC\' %4.4d /", d, y);
          break;
      }
    }

  double 
  dt_tools::ecl2d (const char *buf) const
    {
      char string[512];
      char *ptr_start = 0;
      char *ptr_end = 0;
      int i, year, day, month;
      
      if (!buf)
        return -3;
      
      // Need to add reading check
      if (sscanf (buf, "%d %s %d", &day, string, &year) != 3) 
        {
          return -1;
        }
        
      locale_ucase(string);
      ptr_start = string;
      
      while (*ptr_start == '\'' || *ptr_start == '\"')
        ++ptr_start;
      for (i = 0; ptr_start[i] != '\0'; ++i)
        ;
      ptr_end = &(ptr_start[i]) - 1;
      while (*ptr_end == '\'' || *ptr_end == '\"')
        --ptr_end;
      *(ptr_end + 1) = '\0';    
      if (!strcmp (ptr_start, "JAN"))
        month = 1;
      else if (!strcmp (ptr_start, "FEB"))
        month = 2;
      else if (!strcmp (ptr_start, "MAR"))
        month = 3;
      else if (!strcmp (ptr_start, "APR"))
        month = 4;
      else if (!strcmp (ptr_start, "MAY"))
        month = 5;
      else if (!strcmp (ptr_start, "JUN"))
        month = 6;
      else if (!strcmp (ptr_start, "JLY") || !strcmp (ptr_start, "JUL"))
        month = 7;
      else if (!strcmp (ptr_start, "AUG"))
        month = 8;
      else if (!strcmp (ptr_start, "SEP"))
        month = 9;
      else if (!strcmp (ptr_start, "OCT"))
        month = 10;
      else if (!strcmp (ptr_start, "NOV"))
        month = 11;
      else if (!strcmp (ptr_start, "DEC"))
        month = 12;
      else
        return -4;

      return ymd2d (year, month, day);
    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string
  dt_tools::py_str () const
    {
      std::stringstream s;
      //s << file_name << "\n";
      return s.str ();
    }

  boost::python::list 
  dt_tools::d2date (const double d) const
    {
      boost::python::list l;
      int day, month, year, hour, minute, second;
      double dd;
      dd = d2ymd (d, year, month, day);
      t2hms (dd, hour, minute, second);
      l.append (year);
      l.append (month);
      l.append (day);
      l.append (hour);
      l.append (minute);
      l.append (second);
      return l;
    }

  double 
  dt_tools::date2d (int year, int month, int day, int hour, int minute, int second) const
    {
      return ymd2d (year, month, day) + hms2t (hour, minute, second);
    }

  std::string 
  dt_tools::d2str (double d) const
    {
      char b[1024];
      d2cstr (d, b);
      return std::string (b);
    }
  std::string 
  dt_tools::t2str (double t) const
    {
      char b[1024];
      t2cstr (t, b);
      return std::string (b);
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (dt_tools);
  BLUE_SKY_TYPE_STD_COPY (dt_tools);

  BLUE_SKY_TYPE_IMPL(dt_tools,  dt_tools_iface, "dt_tools", "Date time tools", "Date time tools");

}  // blue_sky namespace
