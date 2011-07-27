/** 
 * @file gis.cpp
 * @brief WELL GIS storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-26
 */

#include <memory.h>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
//#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include "bs_kernel.h"
#include "gis.h"


using namespace std;
using namespace boost;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{

  gis::gis (bs_type_ctor_param) 
    {
      sp_prop = BS_KERNEL.create_object ("prop");
      sp_table = BS_KERNEL.create_object ("table");
      if (!sp_prop || !sp_table)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      
    }
  gis::gis (const gis& rhs) 
    //    : bs_refcounter ()
    {
      *this = rhs;
    }
  int 
  gis::read_ver_info (sp_prop_iface prop, const string &s)
    {
      using qi::str;
      using qi::phrase_parse;
      using ascii::space;
      string s1, s2;
      qi::parse (s.begin (), s.end(), str >> ':' >> str, s1, s2);
      cout << "S1 " << s1 << " S2 "<< s2 <<endl
    }

  int 
  gis::read_from_las_file (const std::string &fname)
    {
      fstream file;
      string s;
      int state = 0;
      const string ver_sec = "~V";
      const string wel_sec = "~W";
      const string par_sec = "~P";
      const string cur_sec = "~C";
      const string oth_sec = "~O";
      const string asc_sec = "~A";

      try 
        {
          file.open (fname.c_str ());
        }
      catch (...)
        {
          //TODO error
          ;
        }
      for (;!file.eof ();)
        {
          getline (file, s);
          trim (s);
          if (s[0] == '#')
            continue;
          else if (s[0] == '~')
            {
              if (s[1] == 'V') // version info
                {
                  state = 1;
                }
              else if (s[1] == 'W') // well info
                {
                  state = 2;
                }
              else if (s[1] == 'P') // well info
                {
                  state = 3;
                }
              else if (s[1] == 'C') // well info
                {
                  state = 4;
                }
              else if (s[1] == 'O') // well info
                {
                  state = 5;
                }
              else if (s[1] == 'A') // well info
                {
                  state = 6;
                }
              else
                {
                  //TODO unknown section
                  cout << s << endl;
                }
              continue;
            }
          else
            {
              if (state == 1)
                read_ver_info (sp_prop, s);
            }


        }
      return 0;
    }

#ifdef BSPY_EXPORTING_PLUGIN
  string 
  gis::py_str () const
    {
      stringstream s;
      t_long i, j;
      s << sp_prop->py_str () << "\n";
      s << sp_table->py_str () << endl;
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (gis);
  BLUE_SKY_TYPE_STD_COPY (gis);

  BLUE_SKY_TYPE_IMPL(gis,  gis_iface, "gis", "GIS storage", "realization of well GIS storage");

}  // blue_sky namespace
