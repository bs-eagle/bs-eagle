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
  
  class BS_API_PLUGIN prop : public prop_iface
    {
      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      typedef std::list<std::string>                  list_t;
      typedef boost::archive::binary_iarchive           tia_t;
      typedef boost::archive::binary_oarchive           toa_t;

      // destructor
      virtual ~prop ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:
      //! add new floating point property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_f (const t_double def_value, const std::string &short_name, const std::string &description)
        {
          fp_impl.add (def_value, short_name, description);
        }
      
      //! add new integer property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_i (const t_long def_value, const std::string &short_name, const std::string &description) 
        {
          i_impl.add (def_value, short_name, description);
        }

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_s (const std::string def_value, const std::string &short_name, const std::string &description) 
        {
          s_impl.add (def_value, short_name, description);
        }

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_b (const bool def_value, const std::string &short_name, const std::string &description) 
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
      virtual t_double get_f (const std::string &name) const
        {
          return fp_impl.get (name);
        }

      //! return property value
      virtual t_long get_i (const std::string &name) const
        {
          return i_impl.get (name);
        }

      //! return property value
      virtual std::string get_s (const std::string &name) const
        {
          return s_impl.get (name);
        }

      //! return property value
      virtual bool get_b (const std::string &name) const
        {
          return b_impl.get (name);
        }

      //! set value
      virtual void set_f (const std::string &name, const t_double value)
        {
          fp_impl.set (name, value);
        }

      //! set value
      virtual void set_i (const std::string &name, const t_long value)
        {
          i_impl.set (name, value);
        }

      //! set value
      virtual void set_s (const std::string &name, const std::string value)
        {
          s_impl.set (name, value);
        }

      //! set value
      virtual void set_b (const std::string &name, const bool value)
        {
          b_impl.set (name, value);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_f (const std::string &name) const
        {
          return fp_impl.check (name);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_i (const std::string &name) const 
        {
          return i_impl.check (name);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_s (const std::string &name) const
        {
          return s_impl.check (name);
        }

      //! check (return false if property set by default, true otherwise
      virtual bool check_b (const std::string &name) const
        {
          return b_impl.check (name);
        }

      //! reset to default value
      virtual void reset_f (const std::string &name)
        {
          fp_impl.reset (name);
        }

      //! reset to default value
      virtual void reset_i (const std::string &name)
        {
          i_impl.reset (name);
        }

      //! reset to default value
      virtual void reset_s (const std::string &name)
        {
          s_impl.reset (name);
        }

      //! reset to default value
      virtual void reset_b (const std::string &name)
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

      virtual void save (toa_t &ar) const
        {
          ar & fp_impl;
          ar & i_impl;
          ar & s_impl;
          ar & b_impl;
        }
      virtual void load (tia_t &ar)
        {
          ar & fp_impl;
          ar & i_impl;
          ar & s_impl;
          ar & b_impl;
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

      list_t get_names_i () const
        {return i_impl.get_names ();}
      list_t get_names_f () const
        {return fp_impl.get_names ();}
      list_t get_names_b () const
        {return b_impl.get_names ();}
      list_t get_names_s () const
        {return s_impl.get_names ();}

      t_double get_def_val_f (const std::string &name) const
        {return fp_impl.get_def_val (name);}
      t_long get_def_val_i (const std::string &name) const
        {return i_impl.get_def_val (name);}
      bool get_def_val_b (const std::string &name) const
        {return b_impl.get_def_val (name);}
      std::string get_def_val_s (const std::string &name) const
        {return s_impl.get_def_val (name);}

      std::string get_description_f (const std::string &name) const
        {return fp_impl.get_description (name);}
      std::string get_description_i (const std::string &name) const
        {return i_impl.get_description (name);}
      std::string get_description_b (const std::string &name) const
        {return b_impl.get_description (name);}
      std::string get_description_s (const std::string &name) const
        {return s_impl.get_description (name);}


#endif //BSPY_EXPORTING_PLUGIN
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      prop_impl <t_double> fp_impl;
      prop_impl <t_long>  i_impl;
      prop_impl <std::string>  s_impl;
      prop_impl <bool>  b_impl;

      BLUE_SKY_TYPE_DECL (prop);
    };

}; //end of blue_sky namespace
#endif //__PROP_H

