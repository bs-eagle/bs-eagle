/**
 * @file bos_reader_iface.h
 * @brief interface for the BOS ascii file reader 
 * @author Oleg Borschuk
 * @version
 * @date 2012-03-01
 */
#ifndef BOS_READER_IFACE_N7NU84AM

#define BOS_READER_IFACE_N7NU84AM

#include <string>

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python/list.hpp>
#endif //BSPY_EXPORTING_PLUGIN

#include "bs_object_base.h"
#include "conf.h"
#include "dt_tools_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()

#define CHAR_BUF_LEN 4096
#define FREAD_CONVERT_CASE 1
#define FREAD_DONT_CONVERT_CASE 0

namespace blue_sky
{
class bos_reader_iface : public objbase
  {
    public:
      typedef smart_ptr<dt_tools_iface, true>                     sp_dt_tools;
    //////////////////////////////////////////
    // METHODS
    // ---------------------------------------

    public:
      /**
       * @brief destructor
       */
      virtual ~bos_reader_iface ()
        {}

      /** 
       * @brief open file for reading 
       * 
       * @param fname   -- <INPUT> file name
       * @param path    -- <INPUT> base path for include files
       * 
       * @return 0 if success, < 0 if error occur
       */
      virtual int open (const char *fname, const char *path) = 0;

      /** 
       * @brief close all opened files
       */
      virtual void close () = 0;

      /** 
       * @brief return string with current position in file
       */
      virtual std::string get_prefix () = 0;

      /** 
       * @brief skip #count floating point values in string
       * 
       * @param start_ptr   -- <INPUT> string pointer
       * @param end_ptr     -- <OUTPUT> string pointer after value skiping 
       * @param count       -- <INPUT> number of skiping values
       * 
       * @return 0 if success
       */
      virtual int skip_fp (char *start_ptr, char **end_ptr, const t_int count) const = 0;

      /** 
       * @brief skip #count integer values in string
       * 
       * @param start_ptr   -- <INPUT> string pointer
       * @param end_ptr     -- <OUTPUT> string pointer after value skiping 
       * @param count       -- <INPUT> number of skiping values
       * 
       * @return 0 if success
       */
      virtual int skip_int (char *start_ptr, char **end_ptr, const t_int count) const = 0;
      
      /** 
       * @brief skip #count words in string
       * 
       * @param start_ptr   -- <INPUT> string pointer
       * @param end_ptr     -- <OUTPUT> string pointer after value skiping 
       * @param count       -- <INPUT> number of skiping values
       * 
       * @return 0 if success
       */
      virtual int skip_str (char *start_ptr, char **end_ptr, const t_int count) const = 0;

      // read line from file
      /** 
       * @brief read line from the opened file, skip comments (starting with '#', ';', '--') 
       * 
       * @param line        -- <OUTPUT> buffer for line
       * @param max_len     -- <INPUT> buffer length
       * @param flg         -- <INPUT> case convertion flag (FREAD_CONVERT_CASE or FREAD_DONT_CONVERT_CASE)
       * 
       * @return < 0 if error occur 
       */
      virtual int read_line (char *line, const int max_len, const int flg = FREAD_CONVERT_CASE) = 0;

      /** 
       * @brief read text block 
       * 
       * @param line        -- <OUTPUT> buffer for text block
       * @param max_len     -- <INPUT> buffer length
       * 
       * @return < 0 if error occur
       */
      virtual int read_text_block (char *line, const int max_len) = 0;

      // For keyword KEY read array from file stream FP to buffer ARRAY.

      /** 
       * @brief read integer array 
       * 
       * @param key         -- <INPUT> keyword for log printing only
       * @param array       -- <OUTPUT> pointer to the array
       * @param len_array   -- <INPUT> array length
       * 
       * @return number of read elements
       */
      virtual t_long read_int_array (const char *key, t_int *array, const t_long len_array) = 0;

      
      /** 
       * @brief read floating point array
       * 
       * @param key         -- <INPUT> keyword for log printing only
       * @param array       -- <OUTPUT> pointer to the array
       * @param len_array   -- <INPUT> array length
       * @param first_flag  -- <INPUT> skip first occurance of '/'
       * 
       * @return number of read elemens  
       */
      virtual t_long read_fp_array (const char *key, t_float *array, const t_long len_array, 
                                    const int first_flag = 0) = 0;

      /** 
       * @brief read floating point table
       * 
       * @param key         -- <INPUT> keyword for log printing only
       * @param dbuf        -- <OUTPUT> allocated buffer for table
       * @param max_len     -- <INPUT> length of the buffer
       * @param num_col     -- <INPUT> number of columns in table
       * 
       * @return number of read values 
       */
      virtual t_long read_fp_table (const char *key, t_float *dbuf,  const t_long max_len, const t_long num_col) = 0;

      /** 
       * @brief unpack string line (5* -> * * * * *)
       * 
       * @param s       -- <INPUT/OUTPUT> string
       * @param flag    -- <INPUT> number of fields to skip 
       * 
       * @return 0 if success
       */
      virtual int unwrap (char *s, const int flag = 0) const = 0;
      
      /** 
       * @brief scan word from string 
       * 
       * @param start_ptr       -- <INPUT> pointer to the start of the string
       * @param end_ptr         -- <OUTPUT> pointer to the rest of the string
       * @param dest            -- <OUTPUT> word
       * @param no_default_p    -- 
       * @param is_default      -- 
       * 
       * @return 0 if success 
       */
      virtual int scanf_s (char *start_ptr, char **end_ptr, char *dest, 
                           int no_default_p = 0, int *is_default = 0) const = 0;
      
      /** 
       * @brief scan word in commas from string
       * 
       * @param start_ptr       -- <INPUT> pointer to the start of the string
       * @param end_ptr         -- <OUTPUT> pointer to the rest of the string
       * @param dest            -- <OUTPUT> word
       * 
       * @return 0 if success 
       */
      virtual int scanf_text (char *start_ptr, char **end_ptr, char *dest) const = 0;
      
      /** 
       * @brief scan integer from string
       * 
       * @param start_ptr       -- <INPUT> pointer to the start of the string
       * @param end_ptr         -- <OUTPUT> pointer to the rest of the string
       * @param dest            -- <OUTPUT> integer
       * @param is_default      -- 
       * 
       * @return 0 if success
       */
      virtual int scanf_int (char *start_ptr, char **end_ptr, t_long *dest, int *is_default = 0) const = 0;
      
      /** 
       * @brief scan floating point from string
       * 
       * @param start_ptr       -- <INPUT> pointer to the start of the string
       * @param end_ptr         -- <OUTPUT> pointer to the rest of the string
       * @param dest            -- <OUTPUT> floating point
       * @param no_default_p    --
       * @param is_default      --
       * 
       * @return 
       */
      virtual int scanf_fp (char *start_ptr, char **end_ptr, t_double *dest, 
                            int no_default_p = 0, int *is_default = 0) const = 0;
      
      /** 
       * @brief scan filename from string
       * 
       * @param dest            -- <OUTPUT> filename
       * @param source          -- <INPUT> string
       * 
       * @return 0 if success
       */
      virtual int scanf_file_name (char *dest, char *source) const = 0;

      /** 
       * @brief skip phrase
       * 
       * @param next        -- <INPUT/OUTPUT> pointer to the string
       * @param delim       -- <INPUT> delimater
       * 
       * @return 
       */
      virtual int get_phrase (char **next, char delim = ';') const = 0;

      /** 
       * @brief read phrase before delimater
       * 
       * @param next        -- <INPUT/OUTPUT> pointer to the string
       * @param buf         -- <OUTPUT> buffer for phrase
       * @param delim       -- <INPUT> delimater
       * 
       * @return 0 if success
       */
      virtual int get_phrase_str (char **next, char *buf, char delim = ';') const = 0;


      /** 
       * @brief read filepath from string before delimater
       * 
       * @param next        -- <INPUT/OUTPUT> pointer to the string
       * @param buf         -- <OUTPUT> buffer for phrase
       * @param delim       -- <INPUT> delimater
       * 
       * @return 0 if success 
       */
      virtual int get_phrase_filepath (char **next, char *buf, char delim = ';') const = 0;


      /** 
       * @brief read integer value from string before delimater
       * 
       * @param next        -- <INPUT/OUTPUT> pointer to the string
       * @param i           -- <OUTPUT> integer
       * @param delim       -- <INPUT> delimater
       * 
       * @return 0 if success  
       */
      virtual int get_phrase_int (char **next, t_long *i, char delim = ';') const = 0;

      /** 
       * @brief read floating point from string before delimater
       * 
       * @param next        -- <INPUT/OUTPUT> pointer to the string
       * @param d           -- <OUTPUT> floating point
       * @param delim       -- <INPUT> delimater
       * 
       * @return 0 if success   
       */
      virtual int get_phrase_double (char **next, t_double *d, char delim = ';') const = 0;

      /** 
       * @brief return dt_tools class
       */
      virtual sp_dt_tools get_dt () = 0;

      /** 
       * @brief move pointer to the next non blanks
       * 
       * @param start_ptr -- <INPUT/OUTPUT> pointer
       */
      virtual void trim_left (char **start_ptr) const = 0;


      /** 
       * @brief convert C string to upper case
       * 
       * @param s   -- <INPUT/OUTPUT> C string
       */
      virtual void locale_ucase (char *s) const = 0;

      /** 
       * @brief return include dir
       */
      virtual std::string get_incdir () const = 0;
#if 0
      /** 
       * @brief read date from string in format 12 'MAR' 2008 
       * 
       * @param buf     -- <INPUT> string 
       * @param date    -- <OUTPUT> date in days
       * 
       * @return 0 if success
       */
      virtual int key_read_date (const char *buf, double &date) const = 0;
#endif //0

#ifdef BSPY_EXPORTING_PLUGIN
      
      /**
       * @brief python print wrapper
       *
       * @return return table description
       */
      virtual std::string py_str () const = 0;

      virtual int init (const std::string &fname, const std::string &path) = 0;
      virtual std::string read_line_str () = 0;

      /** 
       * @brief read lines and check first simbol if it letter 
       * 
       * @return string with keyword
       */
      virtual std::string get_next_keyword () = 0;
#endif //BSPY_EXPORTING_PLUGIN
};

}  // end of bluesky name space


#endif /* end of include guard: BOS_READER_IFACE_N7NU84AM */

