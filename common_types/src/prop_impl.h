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

  // ----------------------
  // METHODS
  // ----------------------
  public:
    // constructor
    prop_impl () {};
    // destructor
    ~prop_impl () {};

    // add new property
    int add (const type_t def_value, const std::string &short_name, const std::string &description);

    //! clear all
    void clear ()
      {
        data.clear ();
      }

    // return property index by name
    int get_index (const std::string &short_name) const;

    //! return property value
    const type_t get (const int idx) const
      {
        return data[idx].value;
      }

    //! set value
    void set (const int idx, const type_t value)
      {
        data[idx].value = value;
        data[idx].flag = true;
      }

    //! check (return false if property set by default, true otherwise
    bool check (const int idx) const
      {
        return data[idx].flag;
      }

    //! reset to default value
    void reset (const int idx)
      {
        data[idx].value = data[idx].def_value;
        data[idx].flag = false;
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
    std::vector<prop_storage> data;
};

#ifdef BSPY_EXPORTING_PLUGIN
template <class type_t> std::string 
prop_impl<type_t>::py_str () const 
{
  std::stringstream s;
  s << "+----------------------------------------------------------------------+\n";
  s << "| Index |    Value   |    Short name    |         Description          |\n";
  for (int i = 0, n = data.size (); i < n; ++i)
    {
          s << "|";
          s.width (7);
          s << i; 
          s << "|";
          s.width (12);
          s << get (i);
          s << "|";
          s.width (18);
          s << data[i].short_name;
          s << "|";
          s.width (30);
          s << data[i].description << "|\n";
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
template <class type_t> int 
prop_impl<type_t>::add (const type_t def_value, const std::string &short_name, const std::string &description)
{
  int i;

  data.push_back (prop_storage ());
  i = data.size () - 1;
  data[i].def_value = data[i].value = def_value;
  data[i].flag = false;
  data[i].short_name = short_name;
  data[i].description = description;

  return i;
}

/** 
 * @brief return property index by given name 
 * 
 * @param short_name    -- <INPUT> given name
 * 
 * @return < 0 if can not find given name
 */
template <class type_t> int 
prop_impl<type_t>::get_index (const std::string &short_name) const
{
  int n = data.size ();
  for (int i = 0; i < n; ++i)
    {
      if (data[i].short_name == short_name)
        return i;
    }
  return -1;
}

/** 
 * @brief reset all properties to default values
 */
template <class type_t> void
prop_impl<type_t>::reset_all ()
{
  int n = data.size ();
  for (int i = 0; i < n; ++i)
    {
      reset (i);
    }
}

#endif //__PROP_IMPL_H

