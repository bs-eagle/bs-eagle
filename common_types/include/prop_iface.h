/** 
 * @file prop_iface.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-27
 */
#ifndef __PROP_IFACE_H
#define __PROP_IFACE_H
#include <string>

#include "bs_object_base.h"
#include "conf.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()


namespace blue_sky
{

  /**
  * \brief properties
  */
  class prop_iface : public objbase
    {
    public:
      typedef std::vector<std::wstring> list_t;

      virtual ~prop_iface ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:
      //! add new floating point property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_f (const t_double def_value, const std::wstring  &short_name, const std::wstring  &description) = 0;
      
      //! add new integer property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_i (const t_long def_value, const std::wstring  &short_name, const std::wstring  &description) = 0; 

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_s (const std::wstring def_value, const std::wstring  &short_name, const std::wstring  &description) = 0; 

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_b (const bool def_value, const std::wstring  &short_name, const std::wstring  &description) = 0; 
 
      //! clear all
      virtual void clear () = 0;

      //! return property value
      virtual t_double get_f (const std::wstring &name) const = 0;

      //! return property value
      virtual t_long get_i (const std::wstring &name) const = 0;

      //! return property value
      virtual std::wstring get_s (const std::wstring &name) const = 0;

      //! return property value
      virtual bool get_b (const std::wstring &name) const = 0;

      //! set value
      virtual void set_f (const std::wstring &name, const t_double  value) = 0;

      //! set value
      virtual void set_i (const std::wstring &name, const t_long  value) = 0;

      //! set value
      virtual void set_s (const std::wstring &name, const std::wstring  value) = 0;

      //! set value
      virtual void set_b (const std::wstring &name, const bool  value) = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_f (const std::wstring &name) const = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_i (const std::wstring &name) const = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_s (const std::wstring &name) const = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_b (const std::wstring &name) const = 0;

      //! reset to default value
      virtual void reset_f (const std::wstring &name) = 0;

      //! reset to default value
      virtual void reset_i (const std::wstring &name) = 0;

      //! reset to default value
      virtual void reset_s (const std::wstring &name) = 0;

      //! reset to default value
      virtual void reset_b (const std::wstring &name) = 0;

      //! reset all
      virtual void reset_all () = 0;

      /** 
       * @brief pack(serialize) all information of class to text string 
       * 
       * @return string
       */
      virtual std::string to_str () const = 0; 

      /** 
       * @brief Reastore all class information from input text string
       * 
       * @param s -- <INPUT> string
       */
      virtual void from_str (const std::string &s) = 0;

      virtual list_t get_names_i () const = 0;
      virtual list_t get_names_f () const = 0;
      virtual list_t get_names_b () const = 0;
      virtual list_t get_names_s () const = 0;

      virtual t_double get_def_val_f (const std::wstring &name) const = 0;
      virtual t_long get_def_val_i (const std::wstring &name) const = 0;
      virtual bool get_def_val_b (const std::wstring &name) const = 0;
      virtual std::wstring get_def_val_s (const std::wstring &name) const = 0;

      virtual std::wstring get_description_f (const std::wstring &name) const = 0;
      virtual std::wstring get_description_i (const std::wstring &name) const = 0;
      virtual std::wstring get_description_b (const std::wstring &name) const = 0;
      virtual std::wstring get_description_s (const std::wstring &name) const = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::wstring py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

    };
};

#endif //__PROP_IFACE_H


