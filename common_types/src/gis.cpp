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

  namespace client
  {
    namespace qi = boost::spirit::qi;

    typedef std::pair<std::string, std::string> pairs_type;

    template <typename Iterator>
    struct key_value_sequence_ordered 
      : qi::grammar<Iterator, pairs_type()>
    {
        key_value_sequence_ordered()
          : key_value_sequence_ordered::base_type(pair)
        {
            pair  =  key >> -(':' >> value);
            key   =  *(qi::char_ - qi::char_(":"));
            value = *qi::char_;
        }

        qi::rule<Iterator, std::pair<std::string, std::string>()> pair;
        qi::rule<Iterator, std::string()> key, value;
    };

    template <typename Iterator>
    struct key_value_sequence_ordered_dot 
      : qi::grammar<Iterator, pairs_type()>
    {
        key_value_sequence_ordered_dot()
          : key_value_sequence_ordered_dot::base_type(pair)
        {
            pair  =  key >> -('.' >> value);
            key   =  *(qi::char_ - qi::char_("."));
            value = *qi::char_;
        }

        qi::rule<Iterator, std::pair<std::string, std::string>()> pair;
        qi::rule<Iterator, std::string()> key, value;
    };
  }



  int 
  gis::read_ver_info (sp_prop_iface prop, std::string &s)
    {
      namespace qi = boost::spirit::qi;

      std::string::iterator begin = s.begin();
      std::string::iterator end = s.end();

      client::key_value_sequence_ordered<std::string::iterator> p;
      client::key_value_sequence_ordered_dot<std::string::iterator> p_dot;
      client::pairs_type v;
      client::pairs_type v_dot;

      if (!qi::parse(begin, end, p, v))
        {
          std::cout << "Error: Parsing failed\n";
        }
      else
        {
          //std::cout << v.first << "-----" <<v.second << std::endl;
          std::string::iterator b_dot = v.first.begin();
          std::string::iterator e_dot = v.first.end();
          if (!qi::parse(b_dot, e_dot, p_dot, v_dot))
            {
              std::cout << "Error: Parsing failed\n";
            }
          else
            {
              trim (v_dot.first);
              trim (v_dot.second);
              //std::cout << v_dot.first << "-----" <<v_dot.second << std::endl;
              prop->add_property_s ("", v_dot.first, v.second);
              prop->set_s (v_dot.first, v_dot.second);
            }
           
        }
      return 0;

    }
  float str2float (std::string &s)
    {
      namespace qi = boost::spirit::qi;
      float f;
      std::string::iterator b_f = s.begin();
      std::string::iterator e_f = s.end();

      if (!qi::parse (b_f, e_f, qi::float_, f))
        {
          std::cout << "Error: Parsing failed\n";
          return 0;
        }
      return f;
    }

  int 
  gis::read_wel_info (sp_prop_iface prop, std::string &s)
    {
      namespace qi = boost::spirit::qi;

      std::string::iterator begin = s.begin();
      std::string::iterator end = s.end();

      client::key_value_sequence_ordered<std::string::iterator> p;
      client::key_value_sequence_ordered_dot<std::string::iterator> p_dot;
      client::pairs_type v;
      client::pairs_type v_dot;

      if (!qi::parse(begin, end, p, v))
        {
          std::cout << "Error: Parsing failed\n";
        }
      else
        {
          //std::cout << v.first << "-----" <<v.second << std::endl;
          std::string::iterator b_dot = v.first.begin();
          std::string::iterator e_dot = v.first.end();
          if (!qi::parse(b_dot, e_dot, p_dot, v_dot))
            {
              std::cout << "Error: Parsing failed\n";
            }
          else
            {
              trim (v_dot.first);
              trim (v_dot.second);
              trim (v.second);
              if (v_dot.second == "")
                {
                  prop->add_property_s ("", v_dot.first, "");
                  prop->set_s (v_dot.first, v.second);
                }
              else
                {
                  //std::cout << v_dot.first << "-----" <<v_dot.second << std::endl;
                  prop->add_property_f (0, v_dot.first, v.second);
                  prop->set_f (v_dot.first, str2float (v_dot.second));
                }
            }
           
        }
      return 0;

    }
  int 
  gis::read_par_info (sp_prop_iface prop, std::string &s)
    {
      namespace qi = boost::spirit::qi;

      std::string::iterator begin = s.begin();
      std::string::iterator end = s.end();

      client::key_value_sequence_ordered<std::string::iterator> p;
      client::key_value_sequence_ordered_dot<std::string::iterator> p_dot;
      client::pairs_type v;
      client::pairs_type v_dot;

      if (!qi::parse(begin, end, p, v))
        {
          std::cout << "Error: Parsing failed\n";
        }
      else
        {
          //trim (v.first);
          //trim (v.second);

          std::string::iterator b_dot = v.first.begin();
          std::string::iterator e_dot = v.first.end();

          if (!qi::parse(b_dot, e_dot, p_dot, v_dot))
            {
              std::cout << "Error: Parsing failed\n";
            }
          else
            {
              trim (v_dot.first);
              trim (v_dot.second);
              if (v_dot.second == "")
                {
                  prop->add_property_s ("", v.first, "");
                  prop->set_s (v.first, v.second);
                }
              else
                {
                  prop->add_property_f (0, v.first, "");
                  prop->set_f (v.first, str2float (v.second));

                }

            }

        }
      return 0;
    }
  int 
  gis::read_cur_info (sp_prop_iface prop, std::string &s, int n)
    {
      namespace qi = boost::spirit::qi;

      std::string::iterator begin = s.begin();
      std::string::iterator end = s.end();

      client::key_value_sequence_ordered<std::string::iterator> p;
      client::key_value_sequence_ordered_dot<std::string::iterator> p_dot;
      client::pairs_type v;
      client::pairs_type v_dot;

      if (!qi::parse(begin, end, p, v))
        {
          std::cout << "Error: Parsing failed\n";
        }
      else
        {
          trim (v.first);
          trim (v.second);
          std::string name = std::string ("param") 
                             + boost::lexical_cast<std::string> (n);
          prop->add_property_s ("", name, v.second);
          prop->set_s (name, v.first);
        }
      return 0;
    }
  int 
  gis::read_asc_info (std::vector<float> &v, std::string &s)
    {
      namespace qi = boost::spirit::qi;
      using boost::spirit::ascii::space;

      std::string::iterator begin = s.begin();
      std::string::iterator end = s.end();

      if (!qi::phrase_parse (begin, end, *qi::double_, space, v))
        {
          std::cout << "Error: Parsing failed\n";
        }
      //std::cout << v.size () << std::endl;
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
                  double stop, start, step;
                  int n;
                  start = sp_prop->get_f ("STRT");
                  stop = sp_prop->get_f ("STOP");
                  step = sp_prop->get_f ("STEP"); 
                  n = (int)((stop - start) / step);
                  sp_table->init (n, param_counter);
                  for (int i = 0; i < param_counter; ++i)
                    {
                      sp_table
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
                read_wel_info (sp_prop, s);
              else if (state == 3)
                read_par_info (sp_prop, s);
              else if (state == 4)
                {
                  read_cur_info (sp_prop, s, param_counter);
                  ++param_counter;
                }
              else if (state == 6)
                {
                  std::vector <float> v;
                  read_asc_info (v, s);
                }

            }


        }
      return 0;
    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  gis::py_str () const
    {
      std::stringstream s;
      s << sp_prop->py_str () << "\n";
      //s << sp_table->py_str () << std::endl;
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (gis);
  BLUE_SKY_TYPE_STD_COPY (gis);

  BLUE_SKY_TYPE_IMPL(gis,  gis_iface, "gis", "GIS storage", "realization of well GIS storage");

}  // blue_sky namespace
