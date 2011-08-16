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
#include <vector>

#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include "bs_kernel.h"
#include "gis.h"


using namespace boost;

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
        : bs_refcounter ()
    {
      *this = rhs;
    }

  void split_str (std::string &s,
                  std::string &name,
                  std::string &units,
                  std::string &data,
                  std::string &description)
    {
      size_t colon_pos;
      size_t dot_pos; 
      size_t space_pos;

      colon_pos = s.find_last_of (':');
      description = s.substr (colon_pos + 1);
      s.erase (colon_pos);
      dot_pos = s.find ('.');
      if (dot_pos == s.npos) // dot not found
        {
          name = s;
          trim (name);
        }
      else
        {
          units = s.substr (dot_pos + 1);
          s.erase (dot_pos);
          name = s;
          trim (name);
          space_pos = units.find (' ');
          if (space_pos != units.npos)
            {
              data = units.substr (space_pos + 1);
              units.erase (space_pos);
              trim(data);
            }
          trim (units);
        }
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
  gis::read_ver_info (sp_prop_iface prop, std::string &s)
    {
      std::string description;
      std::string name;
      std::string units;
      std::string data;

      split_str (s, name, units, data, description);
      if (name == "VERS")
        {
          prop->add_property_f (0, name, description);
          prop->set_f (name, str2T<double>(data));
        }
      else
        {
          prop->add_property_s ("", name, description);
          prop->set_s (name, data);
        }
      return 0;
    }

  int 
  gis::read_wel_info (sp_prop_iface prop, std::string &s, double ver)
    {
      std::string description;
      std::string name;
      std::string units;
      std::string data;

      split_str (s, name, units, data, description);
      if (name == "STRT"
          || name == "STOP"
          || name == "STEP"
          || name == "NULL")
        {
          prop->add_property_f (0, name, units);
          prop->set_f (name, str2T<double> (data));
        }
      else if (name == "LIC")
        {
          prop->add_property_i (0, name, description);
          prop->set_i (name, str2T<int> (data));
        }
      else
        {
          if (ver <= 1.2 && data == "")
            {
              prop->add_property_s ("", name, "");
              prop->set_s (name, description);
            }
          else if (data != "")
            {
              prop->add_property_s ("", name, description);
              prop->set_s (name, data);
            }
        }
      return 0;
    }

  int 
  gis::read_par_info (sp_prop_iface prop, std::string &s, double ver)
    {
      std::string description;
      std::string name;
      std::string units;
      std::string data;

      split_str (s, name, units, data, description);
      if (ver <= 1.2 && data == "")
        {
          prop->add_property_s ("", name, "");
          prop->set_s (name, description);
        }
      else if (data != "")
        {
          prop->add_property_s ("", name, description);
          prop->set_s (name, data);
        }
      return 0;
    }

  int 
  gis::read_cur_info (sp_prop_iface prop, std::string &s, int n)
    {
      std::string description;
      std::string name;
      std::string units;
      std::string data;

      std::string param = std::string ("param") 
                         + boost::lexical_cast<std::string> (n);
      split_str (s, name, units, data, description);

      if (units != "")
        name = name + " (" + units + ")";
      prop->add_property_s ("", param, description);
      prop->set_s (param, name);
      return 0;
    }
  int 
  gis::read_asc_info (std::vector<t_double> &v, std::string &s, 
                      int n, std::fstream &file)
    {
      namespace qi = boost::spirit::qi;
      using boost::spirit::ascii::space;

      for (;!file.eof ();)
        {
          std::string::iterator begin = s.begin();
          std::string::iterator end = s.end();
          if (!qi::phrase_parse (begin, end, *qi::double_, space, v))
            {
              std::cout << "Error: Parsing failed double\n";
            }
          if ((int)v.size () < n)
            {
              if (sp_prop->get_s ("WRAP") != "YES")
                {
                  std::cerr << "Error: wrong number of parameters " << v.size () 
                            << " should be " << n << std::endl;
                  return -1;
                }
              std::getline (file, s);
            }
          else
            {
              break;
            }
        }
      //std::cout << v.size () << n << std::endl;
      return 0;
    }

  int 
  gis::read_from_las_file (const std::string &fname)
    {
      std::fstream file;
      std::string s;
      int state = 0;
      const std::string ver_sec = "~V";
      const std::string wel_sec = "~W";
      const std::string par_sec = "~P";
      const std::string cur_sec = "~C";
      const std::string oth_sec = "~O";
      const std::string asc_sec = "~A";
      int param_counter = 0;
      double ver = 0;
      std::vector <t_double> v;

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
          std::getline (file, s);
          trim (s);
          if (s[0] == '#' || s == "")
            continue;
          else if (s[0] == '~')
            {
              if (s[1] == 'V') // version info
                {
                  state = 1;
                }
              else if (s[1] == 'W') // well info
                {
                  ver = sp_prop->get_f ("VERS");
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
                  sp_table->init (0, param_counter);
                  for (int i = 0; i < param_counter; ++i)
                    {
                      std::string name = std::string ("param") 
                                         + boost::lexical_cast<std::string> (i);
                      sp_table->set_col_name (i, sp_prop->get_s (name));
                    }
                }
              else
                {
                  //TODO unknown section
                std::cout << s << std::endl;
                }
              continue;
            }
          else
            {
              if (state == 1)
                read_ver_info (sp_prop, s);
              else if (state == 2)
                read_wel_info (sp_prop, s, ver);
              else if (state == 3)
                read_par_info (sp_prop, s, ver);
              else if (state == 4)
                {
                  read_cur_info (sp_prop, s, param_counter);
                  ++param_counter;
                }
              else if (state == 6)
                {
                  v.clear ();
                  int ret_code = read_asc_info (v, s, param_counter, file);
                  if (ret_code)
                    return -1;
                  sp_table->push_back (v);
                }

            }


        }
      return 0;
    }

  gis::sp_gis_t 
  gis::check_serial () const
    {
      std::ostringstream oss;
      std::istringstream iss;

      boost::archive::text_oarchive oar(oss);

      sp_gis_t sp_gis = BS_KERNEL.create_object ("gis");

      save (oar);
      iss.str (oss.str ());
      boost::archive::text_iarchive iar(iss);
      sp_gis->load (iar);
      return sp_gis;

    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  gis::to_str () const
    {
      std::ostringstream oss;

      boost::archive::text_oarchive oar(oss);

      save (oar);
      return oss.str ();
    }
  void 
  gis::from_str (const std::string &s)
    {
      std::istringstream iss;

      iss.str (s);
      boost::archive::text_iarchive iar(iss);
      load (iar);
    }
  std::string 
  gis::py_str () const
    {
      std::stringstream s;
      s << sp_prop->py_str () << "\n";
      s << sp_table->py_str () << std::endl;
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (gis);
  BLUE_SKY_TYPE_STD_COPY (gis);

  BLUE_SKY_TYPE_IMPL(gis,  gis_iface, "gis", "GIS storage", "realization of well GIS storage");

}  // blue_sky namespace
