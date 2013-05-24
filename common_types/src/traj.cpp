/** 
 * @file traj.cpp
 * @brief wellbore trajectory implementation
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-28
 */

#include <memory.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>


#include "bs_kernel.h"
#include "traj.h"
#include "bs_misc.h"

using namespace boost;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{
  traj::traj (bs_type_ctor_param) 
    {
      sp_table = BS_KERNEL.create_object ("table");
      if (!sp_table)
        {
          bs_throw_exception ("Type (table) not registered");
        }
      
    }
  traj::traj (const traj& rhs) 
        : bs_refcounter ()
    {
      *this = rhs;
    }

  template <class T>
  T str2T (std::string &s)
    {
      T f;
      trim (s);
      try
        {
          f = boost::lexical_cast<T>(s);
        }
      catch (boost::bad_lexical_cast &)
        {
          f = 0;
        }
      return f;
    }

  int 
  traj::read_from_dev_file (const std::string &fname)
    {
      std::fstream file;
      std::string s;
      std::vector <double> v;
      int state = 0;
      std::vector<std::string> sv;
      namespace fs = boost::filesystem;

      fs::path p(fname);
      try
        {
          if (!fs::exists (p) || !is_regular_file (p))
            {
# ifdef BOOST_POSIX_API
              std::cout  << "Error: can not open file " << p << std::endl;
# else  // BOOST_WINDOWS_API
              //std::cout << "Error: can not open file " << p.file_string() << std::endl;
# endif
              return -2;
            }
        }
      catch (const fs::filesystem_error& ex)
        {
          std::cout << ex.what() << '\n';
          return -1;
        }

      try 
        {
          std::string sss = p.string ();
          std::cout << "Read DEV file: " << sss << std::endl;
          file.open (sss.c_str ());
        }
      catch (...)
        {
          fprintf (stderr, "Error: cannot open file %s\n", p.string ().c_str ());
          return -2;
        }
      for (;!file.eof ();)
        {
          std::getline (file, s);
          trim (s);
          if (s.length() < 1 || s[0] == '#')
            continue;
          if (state == 0) // read header
            {
              boost::split (sv, s, boost::is_any_of("\t "), boost::token_compress_on);
              sp_table->init (0, sv.size ());
              for (size_t i = 0; i < sv.size (); ++i)
			    {
                  sp_table->set_col_name (i, str2wstr(sv[i]));
			    }
              state = 1;
            }
          else
            {
              namespace qi = boost::spirit::qi;
              using boost::spirit::ascii::space;
              std::string::iterator begin = s.begin();
              std::string::iterator end = s.end();
              v.clear ();
              if (!qi::phrase_parse (begin, end, *qi::double_, space, v))
                {
                  std::cout << "Error: Parsing failed double\n";
                }
              else
                {
                  if (v.size () != sv.size ())
                    continue;
                  sp_table->push_back (v);

                }
              
            }
        }
      return 0;
    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  traj::to_str () const
    {
      std::ostringstream oss;

      boost::archive::text_oarchive oar(oss);

      save (oar);
      return oss.str ();
    }
  void 
  traj::from_str (const std::string &s)
    {
      std::istringstream iss;

      iss.str (s);
      boost::archive::text_iarchive iar(iss);
      load (iar);
    }
  std::string 
  traj::py_str () const
    {
      std::stringstream s;
      s << sp_table->py_str () << std::endl;
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (traj);
  BLUE_SKY_TYPE_STD_COPY (traj);

  BLUE_SKY_TYPE_IMPL(traj,  traj_iface, "traj", "Wellbore trajectory storage", "realization of wellbore trajectory storage");

}
