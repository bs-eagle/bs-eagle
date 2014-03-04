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
#include "bs_serialize_decl.h"

namespace blue_sky
{
  
  class BS_API_PLUGIN prop : public prop_iface
    {
      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      typedef std::list<std::wstring>                  list_t;

      // destructor
      virtual ~prop ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:
      //! add new floating point property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_f (const t_double def_value, const std::wstring &short_name, const std::wstring &description)
        {
          fp_impl.add (def_value, short_name, description);
        }

      //! add new integer property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_i (const t_long def_value, const std::wstring &short_name, const std::wstring &description) 
        {
          i_impl.add (def_value, short_name, description);
        }

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_s (const std::wstring def_value, const std::wstring &short_name, const std::wstring &description) 
        {
          s_impl.add (def_value, short_name, description);
        }

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_b (const bool def_value, const std::wstring &short_name, const std::wstring &description) 
        {
          b_impl.add (def_value, short_name, description);
        }
 
      //! clear all
      virtual void clear ()
        {
          fp_impl.clear ();
          i_impl.clear ();
          s_impl.clear ();
          b_impl.clear ();
        }

      //! return property value
      virtual t_double get_f (const std::wstring &name) const
        {
          return fp_impl.get (name);
        }

      //! return property value
      virtual t_long get_i (const std::wstring &name) const
        {
          return i_impl.get (name);
        }

      //! return property value
      virtual std::wstring get_s (const std::wstring &name) const
        {
          return s_impl.get (name);
        }

      //! return property value
      virtual bool get_b (const std::wstring &name) const
        {
          return b_impl.get (name);
        }

      //! set value
      virtual void set_f (const std::wstring &name, const t_double value)
        {
          fp_impl.set (name, value);
        }

      //! set value
      virtual void set_i (const std::wstring &name, const t_long value)
        {
          i_impl.set (name, value);
        }

      //! set value
      virtual void set_s (const std::wstring &name, const std::wstring value)
        {
          s_impl.set (name, value);
        }

      //! set value
      virtual void set_b (const std::wstring &name, const bool value)
        {
          b_impl.set (name, value);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_f (const std::wstring &name) const
        {
          return fp_impl.check (name);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_i (const std::wstring &name) const 
        {
          return i_impl.check (name);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_s (const std::wstring &name) const
        {
          return s_impl.check (name);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_b (const std::wstring &name) const
        {
          return b_impl.check (name);
        }

      //! reset to default value
      virtual void reset_f (const std::wstring &name)
        {
          fp_impl.reset (name);
        }

      //! reset to default value
      virtual void reset_i (const std::wstring &name)
        {
          i_impl.reset (name);
        }

      //! reset to default value
      virtual void reset_s (const std::wstring &name)
        {
          s_impl.reset (name);
        }

      //! reset to default value
      virtual void reset_b (const std::wstring &name)
        {
          b_impl.reset (name);
        }

      //! reset all
      virtual void reset_all ()
        {
          fp_impl.reset_all ();
          i_impl.reset_all ();
          s_impl.reset_all ();
          b_impl.reset_all ();
        }

      /** 
       * @brief pack(serialize) all information of class to text string 
       * 
       * @return string
       */
      virtual std::string to_str () const; 

      /** 
       * @brief Reastore all class information from input text string
       * 
       * @param s -- <INPUT> string
       */
      virtual void from_str (const std::string &s);

      list_t get_names_i () const
        {return i_impl.get_names ();}
      list_t get_names_f () const
        {return fp_impl.get_names ();}
      list_t get_names_b () const
        {return b_impl.get_names ();}
      list_t get_names_s () const
        {return s_impl.get_names ();}

      t_double get_def_val_f (const std::wstring &name) const
        {return fp_impl.get_def_val (name);}
      t_long get_def_val_i (const std::wstring &name) const
        {return i_impl.get_def_val (name);}
      bool get_def_val_b (const std::wstring &name) const
        {return b_impl.get_def_val (name);}
      std::wstring get_def_val_s (const std::wstring &name) const
        {return s_impl.get_def_val (name);}

      std::wstring get_description_f (const std::wstring &name) const
        {return fp_impl.get_description (name);}
      std::wstring get_description_i (const std::wstring &name) const
        {return i_impl.get_description (name);}
      std::wstring get_description_b (const std::wstring &name) const
        {return b_impl.get_description (name);}
      std::wstring get_description_s (const std::wstring &name) const
        {return s_impl.get_description (name);}

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::wstring py_str () const
        {
          return L"";
        }
#endif //BSPY_EXPORTING_PLUGIN
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      prop_impl <t_double> fp_impl;
      prop_impl <t_long>  i_impl;
      prop_impl <std::wstring>  s_impl;
      prop_impl <bool>  b_impl;

      BLUE_SKY_TYPE_DECL (prop);
      friend class bs_serialize;
    };

}; //end of blue_sky namespace
#endif //__PROP_H

