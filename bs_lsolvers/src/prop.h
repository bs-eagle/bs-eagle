/** 
 * @file prop.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-27
 */
#ifndef __PROP_H
#define __PROP_H
#include <string>
#include <sstream>

#include "prop_iface.h"
#include "prop_impl.h"

namespace blue_sky
{
  
  template <class t_double, class i_type_t, class s_type_t, class b_type_t>
  class BS_API_PLUGIN prop : public prop_iface<t_double, i_type_t, s_type_t, b_type_t>
    {
      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      ~prop ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:
      //! add new floating point property to the list (return index of the property or < 0 if error occur)
      virtual int add_property_f (const t_double def_value, const std::string &short_name, const std::string &description)
        {
          return fp_impl.add (def_value, short_name, description);
        }
      
      //! add new integer property to the list (return index of the property or < 0 if error occur)
      virtual int add_property_i (const i_type_t def_value, const std::string &short_name, const std::string &description) 
        {
          return i_impl.add (def_value, short_name, description);
        }

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual int add_property_s (const s_type_t def_value, const std::string &short_name, const std::string &description) 
        {
          return s_impl.add (def_value, short_name, description);
        }

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual int add_property_b (const b_type_t def_value, const std::string &short_name, const std::string &description) 
        {
          return b_impl.add (def_value, short_name, description);
        }
 
      //! clear all
      virtual void clear ()
        {
          fp_impl.clear ();
          i_impl.clear ();
          s_impl.clear ();
          b_impl.clear ();
        }

      //! return property index by name
      virtual int get_index_f (const std::string &short_name) const
        {
          return fp_impl.get_index (short_name);
        }

      //! return property index by name
      virtual int get_index_i (const std::string &short_name) const
        {
          return i_impl.get_index (short_name);
        }

      //! return property index by name
      virtual int get_index_s (const std::string &short_name) const
        {
          return s_impl.get_index (short_name);
        }

      //! return property index by name
      virtual int get_index_b (const std::string &short_name) const
        {
          return b_impl.get_index (short_name);
        }

      //! return property value
      virtual const t_double get_f (const int idx) const
        {
          return fp_impl.get (idx);
        }

      //! return property value
      virtual const i_type_t get_i (const int idx) const
        {
          return i_impl.get (idx);
        }

      //! return property value
      virtual const s_type_t get_s (const int idx) const
        {
          return s_impl.get (idx);
        }

      //! return property value
      virtual const b_type_t get_b (const int idx) const
        {
          return b_impl.get (idx);
        }

      //! set value
      virtual void set_f (const int idx, const t_double value)
        {
          fp_impl.set (idx, value);
        }

      //! set value
      virtual void set_i (const int idx, const i_type_t value)
        {
          i_impl.set (idx, value);
        }

      //! set value
      virtual void set_s (const int idx, const s_type_t value)
        {
          s_impl.set (idx, value);
        }

      //! set value
      virtual void set_b (const int idx, const b_type_t value)
        {
          b_impl.set (idx, value);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_f (const int idx) const
        {
          return fp_impl.check (idx);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_i (const int idx) const 
        {
          return i_impl.check (idx);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_s (const int idx) const
        {
          return s_impl.check (idx);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_b (const int idx) const
        {
          return b_impl.check (idx);
        }

      //! reset to default value
      virtual void reset_f (const int idx)
        {
          fp_impl.reset (idx);
        }

      //! reset to default value
      virtual void reset_i (const int idx)
        {
          i_impl.reset (idx);
        }

      //! reset to default value
      virtual void reset_s (const int idx)
        {
          s_impl.reset (idx);
        }

      //! reset to default value
      virtual void reset_b (const int idx)
        {
          b_impl.reset (idx);
        }

      //! reset all
      virtual void reset_all ()
        {
          fp_impl.reset_all ();
          i_impl.reset_all ();
          s_impl.reset_all ();
          b_impl.reset_all ();
        }

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;
          
          s << "FP properties:\n";
          s << fp_impl.py_str ();
          s << "INT properties:\n";
          s << i_impl.py_str ();
          s << "STRING properties:\n";
          s << s_impl.py_str ();
          s << "BOOL properties:\n";
          s << b_impl.py_str ();
          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      prop_impl <t_double> fp_impl;
      prop_impl <i_type_t>  i_impl;
      prop_impl <s_type_t>  s_impl;
      prop_impl <b_type_t>  b_impl;

      BLUE_SKY_TYPE_DECL(prop);
    };

}; //end of blue_sky namespace
#endif //__PROP_H

