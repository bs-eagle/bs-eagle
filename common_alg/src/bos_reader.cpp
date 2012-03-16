/** 
 * @file bos_reader.cpp
 * @brief Reader for the BOS ascii files
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2012-03-01
 */

#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "bos_reader.h"
#include "bs_kernel.h"
//#include "localization.h"

using namespace boost;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{

  bos_reader::bos_reader (bs_type_ctor_param)
    {
      fr_file = 0;
      dt_t = BS_KERNEL.create_object ("dt_tools");
    }
  bos_reader::bos_reader (const bos_reader& /*rhs*/)
        : bs_refcounter ()
    {
      fr_file = 0;
      dt_t = BS_KERNEL.create_object ("dt_tools");
    }
  bos_reader::~bos_reader ()
    {
      if (fr_file)
        delete fr_file;
      fr_file = 0;
    }
  int 
  bos_reader::open (const char *fname, const char *path)
    {
      close ();
      fr_file = new FRead (fname, path);
      if (!fr_file)
        {
          // TODO: report error
          return -1;
        }
      return 0;
    }
  
  void 
  bos_reader::close ()
    {
      if (fr_file)
        delete fr_file;
      fr_file = 0;
    }

  std::string 
  bos_reader::get_prefix ()
    {
      if (fr_file)
        {
          return std::string (fr_file->get_prefix ());
        }
      else
        {
          // TODO: report error
          return std::string ("");
        }
    }
  int 
  bos_reader::skip_fp (char *start_ptr, char **end_ptr, const t_int count) const
    {
      if (fr_file)
        {
          return fr_file->skip_d (start_ptr, end_ptr, count);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }

  int 
  bos_reader::skip_int (char *start_ptr, char **end_ptr, const t_int count) const
    {
      if (fr_file)
        {
          return fr_file->skip_u (start_ptr, end_ptr, count);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  int 
  bos_reader::skip_str (char *start_ptr, char **end_ptr, const t_int count) const
    {
      if (fr_file)
        {
          return fr_file->skip_s (start_ptr, end_ptr, count);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  
  int 
  bos_reader::read_line (char *line, const int max_len, const int flg)
    {
      if (fr_file)
        {
          return fr_file->read_line (line, max_len, flg);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  int 
  bos_reader::read_text_block (char *line, const int max_len)
    {
      if (fr_file)
        {
          return fr_file->read_text_block (line, max_len);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }

  t_long 
  bos_reader::read_int_array (const char *key, t_int *array, const t_long len_array)
    {
      if (fr_file)
        {
          t_long i;
          t_long j;
          const int buf_size = 4096;
          char buf[buf_size];

          if (array == 0)               // Check pointer array
            return -1;
          if (key == 0)                 // Check pointer to keyword
            return -2;
          if (len_array == 0)           // Check array length
            return -3;

          for (i = 0; i < len_array;)
            {
              // Read line to character buffer.
              // If error for reading line -  return error
              if ((fr_file->read_line (buf, buf_size)) < 0)
                {
                  return -2; 
                }
              if (buf[0] == '/')        // End of reading array data
                return i;
              j = fr_file->convert_u (array, len_array, i, buf, key);    // Call function to convert string to
              // int array
              if (j < 0)                // Check for error
                {
                  return -3;
                }
              else
                i += j;           // add number of read double to counter
            }
          return i;
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }

  t_long 
  bos_reader::read_fp_array (const char *key, t_float *array, const t_long len_array, const int first_flag)
    {
      if (fr_file)
        {
          t_long i;
          t_long j;
          const int buf_size = 4096;
          char buf[buf_size];

          if (array == 0)               // Check pointer array
            return -1;
          if (key == 0)                 // Check pointer to keyword
            return -2;
          if (len_array == 0)           // Check array length
            return -3;

          for (i = 0; i < len_array;)
            {
              // Read line to chracter buffer.
              // If error for reading line -  return error
              if ((fr_file->read_line (buf, buf_size)) < 0)
                return -2;
              if (buf[0] == '/' && i == 0 && first_flag)        // End of reading array data
                {
                  continue;
                }
              else if (buf[0] == '/')
                {
                  return i;
                }
              j = fr_file->convert_f (array, len_array, i, buf, key);    // Call function to convert string to
              // double array
              if (j < 0)                // Check for error
                return -3;
              else
                i += j;           // add number of read double to counter
            }
          return i;                     // return number of read double
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }

  t_long 
  bos_reader::read_fp_table (const char *key, t_float *dbuf,  const t_long max_len, const t_long num_col)
    {
      t_long i;
      int skip_flag = 1;
      // read double values string till the end of bufer
      for (i = 0; i < max_len / num_col; i++)
        {
          // if number of values not equal to needed number - break
          if (read_fp_array (key, dbuf + i * num_col,
                                    num_col, skip_flag) != num_col)
            break;
          skip_flag = 0;
        }
      // if end of buffer was reached - return error
      if (i > max_len / num_col)
        {
          return -2;
        }
      return i;
    }
  int 
  bos_reader::unwrap (char *s, const int flag) const
    {
      if (fr_file)
        {
          return fr_file->unwrap (s, flag);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  int 
  bos_reader::scanf_s (char *start_ptr, char **end_ptr, char *dest, int no_default_p, int *is_default) const
    {
      if (fr_file)
        {
          return fr_file->scanf_s (start_ptr, end_ptr, dest, no_default_p, is_default);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  int 
  bos_reader::scanf_text (char *start_ptr, char **end_ptr, char *dest) const
    {
      if (fr_file)
        {
          return fr_file->scanf_text (start_ptr, end_ptr, dest);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  int 
  bos_reader::scanf_int (char *start_ptr, char **end_ptr, t_long *dest, int *is_default) const
    {
      if (fr_file)
        {
          return fr_file->scanf_u (start_ptr, end_ptr, dest, is_default);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  int 
  bos_reader::scanf_fp (char *start_ptr, char **end_ptr, t_double *dest, int no_default_p, int *is_default) const
    {
      if (fr_file)
        {
          return fr_file->scanf_d (start_ptr, end_ptr, dest, no_default_p, is_default);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }
  int 
  bos_reader::scanf_file_name (char *dest, char *source) const
    {
      if (fr_file)
        {
          return fr_file->scanf_file_name (dest, source);
        }
      else
        {
          // TODO: report error
          return -1;
        }
    }

  int
  bos_reader::get_phrase (char **next, char delim) const
  {
    if (!(*next))
      return -1;
    char *start = *next;
    while ((**next) != delim && (**next) != '\0')
      {
        ++(*next);
      }
    if (**next != '\0')
      {
        **next = '\0';
        ++(*next);
      }
    trim_right_s (start);
    trim_left (next);
    return 0;
  }

  int
  bos_reader::get_phrase_str (char **next, char *buf, char delim) const
  {
    trim_left (next);
    char *start = *next;
    if (get_phrase (next, delim))
      return -2;
    //printf ("NAMEMMMM: %s\n", start);
    if (*start == '*' || *start == '\0')
      return 0;
    strcpy (buf, start);
    return 0;
  }

  int
  bos_reader::get_phrase_filepath (char **next, char *buf, char delim) const
  {
    int rc = get_phrase_str (next, buf, delim);
    if (rc)
      return rc;

    size_t n = strlen (buf);
#ifdef BOOST_POSIX_API
    for (size_t i = 0; i < n; ++i)
      {
        if (buf[i] == '\\')
          buf[i] = '/';
      }
#else // WINDOWS
    for (size_t i = 0; i < n; ++i)
      {
        if (buf[i] == '/')
          buf[i] = '\\';
      }
#endif //BOOST_POSIX_API
    return 0;
  }

  int
  bos_reader::get_phrase_int (char **next, t_long *i, char delim) const
  {
    trim_left (next);
    char *start = *next;
    long t;
    if (get_phrase (next, delim))
      return -2;
    if (*start == '*' || *start == '\0')
      return 0;
    if (sscanf (start, "%ld", &t) < 1)
      {
        fprintf (stderr, "Error: can not read int from %s\n", start);
        return -1;
      }
    *i = (t_long)t;
    return 0;
  }

  int
  bos_reader::get_phrase_double (char **next, t_double *d, char delim) const
  {
    trim_left (next);
    char *start = *next;
    double t;
    if (get_phrase (next, delim))
      {
        //printf ("get_phrase kkkkkk\n");
        return -2;
      }
    if (*start == '*' || *start == '\0')
      return 0;
    if (sscanf (start, "%lf", &t) < 1)
      {
        fprintf (stderr, "Error: can not read double from %s\n", start);
        return -1;
      }
    *d = (t_double)t;
    return 0;
  }
  void 
  bos_reader::trim_left (char **start_ptr) const
    {
      while ((**start_ptr == '\t' || **start_ptr == ' ' || **start_ptr == ','
              || **start_ptr == ':') && (**start_ptr != '\0'))
        ++(*start_ptr);
    }
  void 
  bos_reader::locale_ucase (char *s) const
    {
      char *ptr = s;
      for (; *ptr != '\0'; ++ptr)
        *ptr = toupper (*ptr);
    }
  std::string 
  bos_reader::get_incdir () const
    {
      if (fr_file)
        {
          return std::string (fr_file->get_incdir ());
        }
      return std::string ();
    }
#if 0
/*! \brief Nonzero if YEAR is a leap year (every 4 years,\n
except every 100th isn't, and every 400th is).  */
#define IS_LEAP_YEAR(YEAR) \
  ((YEAR) % 4 == 0 && ((YEAR) % 100 != 0 || (YEAR) % 400 == 0))

  int 
  bos_reader::key_read_date (const char *buf, double &date) const
    {
      int ret_code;
      int day, month, year;
      int days_in_current_month;
      static const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
      int i;
      
      if (sscanf (buf, "%d.%d.%d", &day, &month, &year) != 3)
        {
          switch (ret_code = key_read_date_ecl (buf, day, month, year))
          {
            case -1:
              {
                fprintf (stderr, "Error: invalid date format: %s\n", buf);
                return -1;
                }
              case -4:
                {
                  fprintf (stderr, "Error: invalid (nonexistent) date: %s\n", buf);
                  return -1;
                }  
            }
        }
     
      // Check month and year
      
      if (year < 1900)
        {
          fprintf (stderr, "Error: years earlier 1900 not supported: %s\n", buf);
          return -1;
        }
      
      if (month < 1 || month > 12)
        {
          fprintf (stderr, "Error: invalid (nonexistent) date: %s\n", buf);
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
          fprintf (stderr, "Error: invalid (nonexistent) date: %s\n", buf);
          return -1;
        }
      
      return 0;
    }

    int
    bos_reader::key_read_date_ecl (const char *buf, int &day, int &month, int &year) const
    {
      char string[512];
      char *ptr_start = 0;
      char *ptr_end = 0;
      int i;
      
      if (!buf)
        return -3;
      
      // Need to add reading check
      if (sscanf (buf, "%d %s %d", &day, string, &year) != 3) 
        return -1;
        
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
      return 0;
    }
#endif //0

#ifdef BSPY_EXPORTING_PLUGIN
  std::string
  bos_reader::py_str () const
    {
      std::stringstream s;
      //s << file_name << "\n";
      return s.str ();
    }

  int 
  bos_reader::init (const std::string &fname, const std::string &path)
    {
      return open (fname.c_str (), path.c_str ());
    }
  std::string 
  bos_reader::read_line_str ()
    {
      char buf[4096];

      if (read_line (buf, 4096) < 0)
        {
          return std::string ();
        }
      printf ("LINE: %s\n", buf);
      return std::string (buf); 

    }

  std::string 
  bos_reader::get_next_keyword ()
    {
      char buf[4096];
      for (;;)
        {
          if (read_line (buf, 4096) < 0)
            {
              return std::string ();
            }
          if (buf[0] >= 'A' && buf[0] <= 'z')
            {
              return std::string (buf); 
            }
        }
      return std::string ();
    }

#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (bos_reader);
  BLUE_SKY_TYPE_STD_COPY (bos_reader);

  BLUE_SKY_TYPE_IMPL(bos_reader,  bos_reader_iface, "bos_reader", "Reader for BOS ASCII files", "Reader for BOS ASCII files");

}  // blue_sky namespace
