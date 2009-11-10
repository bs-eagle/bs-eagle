#ifndef __READ_CLASS
#define __READ_CLASS

/*!
  \file read_class.h
  \brief Read functions for Keyword Input Language -- class #FRead, list of file included -- class #include_files
*/
#include "array_ext.h"
#include "main_def.h"
#include "err_num_def.h"
#include "throw_exception.h"
#include "str_tools.h"
#include "convert_str_to_array.h"

#include <boost/date_time/gregorian/gregorian.hpp>

namespace blue_sky
{
  #define FREAD_CONVERT_CASE 1
  #define FREAD_DONT_CONVERT_CASE 0

#define DIR_SYMBOL_UNIX '/'
#define DIR_SYMBOL_WIN  '\\'

#ifdef UNIX
  #define DIR_SYMBOL DIR_SYMBOL_UNIX
#else // !UNIX
  #define DIR_SYMBOL DIR_SYMBOL_WIN
#endif //  UNIX

namespace detail {

  inline const char *
  strrchr (const char *str)
  {
    const char *x = ::strrchr (str, DIR_SYMBOL);
    if (x)
      return x;

    x = ::strrchr (str, DIR_SYMBOL_UNIX);
    if (x)
      return x;

    return 0;
  }

} // namespace detail

  /*!
    \class R_FILE
    \brief class stores information about input keyword files
  */
  class BS_API_PLUGIN R_FILE
  {
  public:

    R_FILE ()
    : nstr (0),
    fp (0)
    {
    }

    ~R_FILE ()
    {
      close ();
    }

    //! return file stream
    FILE *
    get_fp () const         
    {
      BS_ASSERT (fp);
      return fp;
    }
    //! return filename
    const std::string &
    get_name () const       
    {
      return name;
    }

    //! open file
    FILE *
    open (const std::string &Name)
    {
      if (!Name.length ())
        bs_throw_exception ("Invalid argument: name is empty");

      if (fp)
        {
          fclose (fp);
        }

      fp = fopen (Name.c_str (), "rt");
      if (!fp)
        bs_throw_exception (boost::format ("Can't open file %s") % Name);

      name = Name;
      return fp;
    }

    void close ()
    {
      if (fp)
        {
          fclose (fp);
        }

      fp = 0;
      nstr = 0;
    }

    size_t        nstr;                 //! position in the file
  private:
    std::string   name;                 //! file name
    FILE          *fp;                  //! file descriptor
  };


  /*!
    \class include_file_node
    \ingroup KeywordLanguage
    \brief include_file_node -- list node for storing file name
  */
  class /*BS_API_PLUGIN*/ include_file 
  {
  public:

    include_file (const std::string &new_file_name = "")
    {
      set_file_name (new_file_name);
    }

    //! set new file name, throw bs_exception if error occur
    void
    set_file_name (const std::string &new_file_name)
    {
      full_file_name = "";
      file_name = "";
      if (!new_file_name.length ())
        {
          bs_throw_exception ("Invalid argument: new_file_name is empty");
        }

      full_file_name = new_file_name;
      const char *last_name = detail::strrchr (full_file_name.c_str ());
      if (last_name)
        {
          file_name = last_name + 1;
        }
      else
        {
          file_name = full_file_name;
        }
    }

    //! return pointer to the file name with full path
    const std::string &
    get_full_file_name () const
    {
      return full_file_name;
    }

    //! return pointer to the file name only
    const std::string &
    get_file_name () const
    {
      return file_name;
    }
  private:
    std::string full_file_name;         //!< file name with full puth
    std::string file_name;              //!< pointer to place in full_file_name variable. Consist only from file name
  };

  /*!
    \class FRead
    \ingroup KeywordLanguage
    \brief FRead is a set of read function needed to read information from
  *               input file, support include files
  */

  class BS_API_PLUGIN FRead : public objbase
    {
      typedef FRead                           this_t;
      typedef smart_ptr<this_t , true>        sp_this_t;
      
    public:
      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL_T(FRead);

    public:
      enum
      {
        FILE_COUNT = 40
      };

      ~FRead ();

      void init (const std::string &file_name, const std::string &dir_name);

      const std::string &
      get_incdir () const;
      
      void close_all ();

      void push (std::string file_name, bool is_main_file = false);
      void pop ();

      boost::gregorian::date read_date(const std::string & str_date);

      //! read line from file
      int read_line (char *line, int max_len, int flg = FREAD_CONVERT_CASE);

      //! read text block from file
      int read_text_block (char *line, int max_len);

      /*!
        \brief For keyword KEY read array from
               file stream fp to buffer ARRAY.\n
               string format: 2*15.8 3*{12 2*5.6 2*{1.1 1.2}} is equals\n
               15.8 15.8 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2
        \param key        Name of calling keyword
        \param array      string buffer
        \param offset     offset in array
        \param len_array  lenght of buffer
        \param first_flag -- if this flag > 0 skip first line starting with '/'
        \return Number of read values in case of success\n
                throw exception in case of reading error
      */
      template <class array_t>
      size_t 
      read_array (const std::string & key, array_t &array, size_t offset = 0,
        size_t len_array = 0, bool first_flag = false)
      {
        boost::array <char, CHAR_BUF_LEN> buf;

        if (len_array == 0)           // Check array length
          {
            len_array = array.size();
            if (len_array == 0)
              {
                bs_throw_exception("Array is empty");
              }
          }

        size_t i = 0;
        for (i = 0; i < len_array;)
          {
            if ((this->read_line (&buf[0], CHAR_BUF_LEN)) < 0)
              {
                bs_throw_exception("Read line error");
              }
            if (buf[0] == '/' && i == 0 && first_flag)        // End of reading array data
              {
                continue;
              }
            else if (buf[0] == '/')
              {
                break;
              }
            size_t j = data_reader::convert (array, offset, i, buf, key, get_prefix ());    
            i += j;
          }

        return i;                     // return number of read items
      }

      //! For keyword KEY read table from file stream FP to buffer DBUF.
      //! Each string contained NUM_COL double values.
      /*!
        \brief For keyword KEY read table from
               file stream this->fp to buffer DBUF.\n
               Each string contained NUM_COL
               double values.
        \param key      Input string of table
        \param table    Array for output
        \param num_col  Number of column in string of table

        \result Number of read values in case of success\n
                throw exception in case of reading error
      */
      template <typename array_t>
      size_t
      read_table (const std::string & key, array_t &table, int num_col)
      {
        typedef typename array_t::value_type item_t;

        seq_vector <item_t> temp;
        temp.resize(num_col);

        size_t line_count = 0;
        if (read_array (key, temp, 0, 0, 1) == (size_t)num_col)
          {
            table.insert (table.end (), temp.begin (), temp.end ());
            line_count++;

            while (this->read_array (key, temp, 0, 0, 0) == (size_t)num_col)
              {
                table.insert (table.end (), temp.begin (), temp.end ());
                line_count++;
              }
          }

        return line_count;
      }

      std::string 
      get_prefix ();

      int eof ()
      {
        return feof (fp);
      }

      typedef seq_vector <include_file> file_list_t;

      //! set pointer to the external include list
      void 
      set_out_inc_list (const file_list_t &ext_list)
      {
        out_inc_list.assign (ext_list.begin (), ext_list.end ());
      }
      
    private:


      FILE                              *fp;
      std::string                       include_dir_;
      file_list_t                       out_inc_list;
      file_list_t                       inc_list;
      int                               top;
      int                               total;
      std::string                       prefix_;
      boost::array <R_FILE, FILE_COUNT> F;
    };
}//bs
#endif//__READ_CLASS
