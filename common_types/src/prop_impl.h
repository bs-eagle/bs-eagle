/** 
 * @file prop_impl.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-27
 */
#ifndef __PROP_IMPL_H
#define __PROP_IMPL_H

#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <list>

#include "throw_exception.h"

template <class type_t>
class prop_impl
{
  //-----------------------
  // TYPE DECL
  // ----------------------
  public:
    struct prop_storage 
      {
        type_t value;
        type_t def_value;
        bool flag;
        std::string short_name;
        std::string description;
      };

    typedef std::map<std::string, prop_storage> map_t;
    typedef std::list<std::string> list_t;
  // ----------------------
  // METHODS
  // ----------------------
  public:
    // constructor
    prop_impl () {};
    // destructor
    ~prop_impl () {};

    // add new property
    void add (const type_t def_value, const std::string &short_name, const std::string &description);

    //! clear all
    void clear ()
      {
        data.clear ();
      }

    //! return property value
    const type_t get (const std::string &name) const
      {
        typename map_t::const_iterator i;

        i = data.find (name);
        if (i != data.end ())
          {
            return i->second.value;
          }
        using namespace blue_sky;
        bs_throw_exception (boost::format ("no property with name: %s") % name);
      }

    //! set value
    void set (const std::string &name, const type_t &value)
      {
        
        typename map_t::iterator i;

        i = data.find (name);
        if (i != data.end ())
          {
            i->second.value = value;
            i->second.flag = true;
          }
        else
          {
            using namespace blue_sky;
            bs_throw_exception (boost::format ("no property with name: %s") % name);
          }
      }

    //! check (return false if property set by default or not exist, true otherwise
    bool check (const std::string &name) const
      {
        typename map_t::const_iterator i;

        i = data.find (name);
        if (i != data.end ())
          {
            return i->second.flag;
          }
        else
          {
            return false;
          }
      }

    //! reset to default value
    void reset (const std::string  &name)
      {
        typename map_t::iterator i;

        i = data.find (name);
        if (i != data.end ())
          {
            i->second.value = i->second.def_value;
            i->second.flag = false;
          }
        else
          {
            using namespace blue_sky;
            bs_throw_exception (boost::format ("no property with name: %s") % name);
          }
      }
    //! return list of avalible names 
    list_t get_names () const 
      {
        typename map_t::const_iterator i, e;
        list_t l;

        for (i = data.begin (), e = data.end (); i != e; ++i)
          {
            l.push_back (i->first);
          }
        return l;
      }

    //! return default value for given name
    const type_t get_def_val (const std::string &name) const
      {
        typename map_t::const_iterator i;

        i = data.find (name);
        if (i != data.end ())
          {
            return i->second.def_value;
          }
        return type_t ();
      }

    //! return default value for given name
    const std::string get_description (const std::string &name) const
      {
        typename map_t::const_iterator i;

        i = data.find (name);
        if (i != data.end ())
          {
            return i->second.description;
          }
        using namespace blue_sky;
        bs_throw_exception (boost::format ("no property with name: %s") % name);
      }


    // reset all
    void reset_all ();

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const;
#endif //BSPY_EXPORTING_PLUGIN

    // -----------------------------------
    // VARIABLES
    // -----------------------------------
  public:
    map_t data;

};

#ifdef BSPY_EXPORTING_PLUGIN
template <class type_t> std::string 
prop_impl<type_t>::py_str () const 
{
  std::stringstream s;
  typename map_t::const_iterator i, e;
  int c = 0;
  s << "+----------------------------------------------------------------------+\n";
  s << "| Index |    Value   |    Short name    |         Description          |\n";
  for (c = 1, i = data.begin (), e = data.end (); i != e; ++i, ++c)
    {
          s << "|";
          s.width (7);
          s << c; 
          s << "|";
          s.width (12);
          s << i->second.value;
          s << "|";
          s.width (18);
          s << i->second.short_name;
          s << "|";
          s.width (30);
          s << i->second.description << "|\n";
    }
  s << "+----------------------------------------------------------------------+\n";
  return s.str ();
}
#endif //BSPY_EXPORTING_PLUGIN
/** 
 * @brief add new property
 * 
 * @param def_value     -- <INPUT> default value 
 * @param short_name    -- <INPUT> name of the property
 * @param description   -- <INPUT> property description
 * 
 * @return 0 if success
 */
template <class type_t> void 
prop_impl<type_t>::add (const type_t def_value, const std::string &short_name, const std::string &description)
{
  typename map_t::iterator i, e;
  std::pair<std::string, prop_storage> p;
  
  i = data.find (short_name);
  if (i == data.end ())
    {
      p.second.def_value = p.second.value = def_value;
      p.second.flag = false;
      p.second.short_name = short_name;
      p.first = short_name;
      p.second.description = description;
      data.insert (p);
    }
}

/** 
 * @brief reset all properties to default values
 */
template <class type_t> void
prop_impl<type_t>::reset_all ()
{
  typename map_t::iterator i, e;

  for (i = data.begin (), e = data.end (); i != e; ++i)
    {
      reset (i->first);
    }
}

#endif //__PROP_IMPL_H

