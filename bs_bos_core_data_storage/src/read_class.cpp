/*!
  \file read_class.cpp
  \brief Read functions for Keyword Input Language class #FRead

  Class #FRead methods.
*/
#include "bs_bos_core_data_storage_stdafx.h"

#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#ifdef UNIX
#include <unistd.h>
#else
#include <io.h>
#endif
#include <fcntl.h>

#include <boost/spirit.hpp>
#include <boost/spirit/core.hpp>
#include <boost/spirit/iterator/file_iterator.hpp>

#include "read_class.h"

#include "main_def.h"
#include "localization.h"

#include "path_tools.h"

namespace blue_sky
  {

    #define SET_LINE_INFO(disp)                                                                   \
    if (disp < 0)                                                                                 \
      {                                                                                           \
        if (feof (this->fp))                                                                      \
          {                                                                                       \
            this->pop ();                                                                         \
            goto beg;                                                                             \
          }                                                                                       \
        else                                                                                      \
          return YS_CANNOT_READ_FROM_FILE;                                                        \
      }                                                                                           \
    IncrLineNumInFile (F, top, disp);

    template <typename file_list_t>
    inline void IncrLineNumInFile (file_list_t &F, size_t top, size_t val)        //!
    {
      if (top > 0)
        {
          F[top - 1].nstr += val;
        }
    }

  FRead::FRead (bs_type_ctor_param)
  {
    fp      = 0;
    top     = 0;
  }

  FRead::FRead (const FRead& src)                     //! Copy Constructor
  : bs_refcounter (src), objbase (src)
  {
    *this = src;
  }

  /*!
    \brief Constructor -- make a new object of class FRead
  *                       and open file 'fname'
    \param fname
    \param dir
  */
  void FRead::init(const std::string &file_name, const std::string &dir_name)
  {
    fp  = 0;
    top = 0;

    out_inc_list.clear ();

    push (file_name, true);
    include_dir_ = path::dirname (dir_name);
  }

  /*!
    \brief Method return pointer to the string of include path
    \return pointer to the string of include path
  */
  const std::string &
  FRead::get_incdir () const
  {
    return include_dir_;
  }

  /*!
    \brief Destructor -- free all allocated memory and close all open files
  */
  FRead::~FRead ()
  {
  }

  /*!
    \brief close all opened files;
  */
  void
  FRead::close_all ()
  {
    if (fp)
      {
        fclose (fp);
        fp = 0;
      }

    // Close all files
    for (int i = 0; i < top; ++i)
      F[i].close ();
    top = 0;
  }

  /*!
    \brief method try to open file with name 'fname', save this->fp in stack
    \param fname -- pointer to the name of file
    \param main_file_p

    \return if success                                       0
    \return if no name or can not open file                 -1
    \return if file already open                            -3
  */
  void
  FRead::push (std::string file_name, bool is_main_file)
  {
    if (!file_name.length ())
      {
        bs_throw_exception ("File name is empty");
      }

    if (!is_main_file && !path::is_absolute_path (file_name))
      {
        file_name = include_dir_ + file_name;
      }

    for (int i = 0; i < top; ++i)
      {
        if (file_name == F[i].get_name ())
          {
            bs_throw_exception (boost::format ("Error in %s: including file %s is already open...") % get_prefix () % file_name);
          }
      }
    if (!(fp = F[top].open (file_name)))
      {
        bs_throw_exception (boost::format ("Error in %s: can't open included file %s") % get_prefix () % file_name);
      }
    top++;

    inc_list.push_back (file_name);
    if (out_inc_list.size ())
      {
        out_inc_list.push_back (file_name);
      }

    if (inc_list.size ())
      {
        BOSOUT (section::read_data, level::medium) << "Include file name = " << file_name << bs_end;
      }
    else
      {
        BOSOUT (section::read_data, level::medium) << "The following files has been included from :" << file_name << bs_end;
      }
  }

  /*!
    \brief method close this->fp file and get open file from stack
    \return if success                                      0
    \return if no open files in stack                       -1
  */
  void
  FRead::pop ()
  {
    if (top - 1 > 0)
      {
        --top;
        F[top].close ();

        fp = F[top - 1].get_fp ();
        BS_ASSERT (fp);
      }
    else
      bs_throw_exception ("File stack is empty");
  }

  boost::gregorian::date FRead::read_date(const std::string & str_date)
  {
    using namespace boost::spirit;

    tm temp_date = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    //Date parsing
    // d.m.y
    rule <> dot_date_p = int_p[assign_a(temp_date.tm_mday)] >> ch_p(".") >>
                         int_p[assign_a(temp_date.tm_mon)] >> ch_p(".") >>
                         int_p[assign_a(temp_date.tm_year)];
    // d 'M' y | d "M" y | d M y
    rule <> month_p = str_p("JAN")[assign_a(temp_date.tm_mon, 0)] |
                      str_p("FEB")[assign_a(temp_date.tm_mon, 1)] |
                      str_p("MAR")[assign_a(temp_date.tm_mon, 2)] |
                      str_p("APR")[assign_a(temp_date.tm_mon, 3)] |
                      str_p("MAY")[assign_a(temp_date.tm_mon, 4)] |
                      str_p("JUN")[assign_a(temp_date.tm_mon, 5)] |
                      str_p("JUL")[assign_a(temp_date.tm_mon, 6)] |
                      str_p("AUG")[assign_a(temp_date.tm_mon, 7)] |
                      str_p("SEP")[assign_a(temp_date.tm_mon, 8)] |
                      str_p("OCT")[assign_a(temp_date.tm_mon, 9)]|
                      str_p("NOV")[assign_a(temp_date.tm_mon, 10)]|
                      str_p("DEC")[assign_a(temp_date.tm_mon, 11)];

    rule <> q_month_p     = month_p 
                              | ch_p("'") >> month_p >> ch_p("'") 
                              | ch_p("\"") >> month_p >> ch_p("\"");

    rule <> space_date_p  = int_p[assign_a(temp_date.tm_mday)] >> +(space_p) 
                              >> q_month_p  >> +(space_p) 
                              >> int_p[assign_a(temp_date.tm_year)];

    rule <> any_date_p    = dot_date_p | space_date_p;

    if (!parse(str_date.c_str(),any_date_p).full)
      throw bs_exception("Fread::read_date","Date parsing failed");

    temp_date.tm_year -= 1900;
    return boost::gregorian::date_from_tm(temp_date);
  }

  /*!
    \brief Read line from this->fp, pass comments begining from '#' or ';',
           work with include files
    \param buf     Name of read buffer
    \param MaxLen  Maximum number of character to read
    \param flg     flag for case conversion if FREAD_CONVERT_CASE convert string to upper case
                                               FREAD_DONT_CONVERT_CASE do not convert string cases


    \return if success                                      Number of readed characters \n
     if no open files in FRead class                 NO_OPEN_FILE \n
     if 'MaxLen' parameter equals zero               0 \n
     if reading error                                -1 \n
     if end of file and nothing have been read       -2\n
     if can not open include file                    -5
  */
  int
  FRead::read_line (char *buf, int MaxLen, int flg)
  {
    const char *inc = GET_TEXT ("include");
    const char *INC = GET_TEXT ("INCLUDE");
    const int inc_len = 7;	// length of word include
    const int f_len = 512;	// max length of file name
    char *word_cand;		// temporary buffer for reading
    int word_len;
    char *fln;			// buffer for include file files name
    int disp;

    if (!buf)
      bs_throw_exception ("Invalid argument: buf is null");

    if (!fp)
      bs_throw_exception ("Invalid state: file not opened");

    if (MaxLen == 0)		// Check for max number of character is not a zero
      return 0;

    while (1)
      {
beg:
        word_cand = buf;
        disp = get_next_not_empty_line (word_cand, MaxLen, this->fp, &word_len);
        SET_LINE_INFO (disp);

        // if exist INCLUDE keyword
        if ((strstr (word_cand, INC) == word_cand)
            || (strstr (word_cand, inc) == word_cand))
          {
            fln = trim_right_s (trim_left_s (word_cand + inc_len, '/'), '/');
            while (*fln == '\0')
              {
                fln = buf;
                disp =
                  get_next_not_empty_line (fln, f_len, this->fp, &word_len);
                SET_LINE_INFO (disp);
                fln = trim_right_s (trim_left_s (fln, '/'), '/');
              }
            this->push (fln);
            continue;
          }
        // pass comment strings
        else if ((word_cand[0] != '#') && (word_cand[0] != ';') &&
                 !(word_cand[0] == '-' && word_cand[1] == '-'))
          {
            word_cand = trim_right_s (word_cand);
            word_len = (int) strlen(word_cand);
            if (!((word_len == 1) && (*word_cand == '/')))
              word_cand = trim_right_s (word_cand, '/');
            break;
          }
      }
    if (flg == FREAD_CONVERT_CASE)
      locale_ucase (buf);
    return word_len;
  }

  /*!
    \brief Read text block up to "\", pass comments begining from '#' or ';' or '--'
    \param buf     Name of read buffer
    \param MaxLen  Maximum number of character to read
    \return if success                               Number of readed characters \n
     if no open files in FRead class                 NO_OPEN_FILE \n
     if 'MaxLen' parameter equals zero               0 \n
     if reading error                                -1 \n
     if end of file and nothing have been read       -2\n
  */
  int
  FRead::read_text_block (char *buf, int MaxLen)
  {
    char *word_cand;		// temporary buffer for reading
    int word_len, len = 0;
    int disp;
    char *end_block;

    if (!buf)
      bs_throw_exception ("Invalid argument: buf is null");

    if (!fp)
      bs_throw_exception ("Invalid state: file not opened");

    if (MaxLen == 0)		// Check for max number of character is not a zero
      return 0;

    word_cand = buf;
    while (1)
      {
beg:
        disp = get_next_not_empty_line (word_cand, MaxLen, this->fp, &word_len);
        SET_LINE_INFO (disp);

        // pass comment strings
        if ((word_cand[0] != '#') && (word_cand[0] != ';') &&
            !(word_cand[0] == '-' && word_cand[1] == '-'))
          {
            end_block = strchr (word_cand, '/');
            if (end_block != NULL)
              {
                int pos = (int)(end_block - word_cand);
                word_cand[pos] = '\0';
                len += pos;
                break;
              }

            word_cand[word_len] = ' ';
            word_cand = word_cand + word_len + 1;
            len += word_len + 1;
          }
      }

    return len;
  }

  /*!
    \brief build and return prefix
    \return Pointer to prefix string
  */
  std::string
  FRead::get_prefix ()
  {
    prefix_ = "";
    for (int i = 0; i < top; ++i)
      {
        prefix_ += (boost::format ("%s '%s' %s %d ") % GET_TEXT ("File") % F[i].get_name () % GET_TEXT ("Line") % F[i].nstr).str ();
      }

    return prefix_;
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(FRead)
  BLUE_SKY_TYPE_STD_COPY(FRead)
  BLUE_SKY_TYPE_IMPL_SHORT(FRead, objbase, "BOS_Core class FRead is a set of read function needed to read information from input file, support include files")



}//blue_sky
