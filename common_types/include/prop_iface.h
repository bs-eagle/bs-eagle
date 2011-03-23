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

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()


namespace blue_sky
{

  /**
  * \brief properties
  */
  template <class fp_type_t, class i_type_t, class s_type_t, class b_type_t>
  class prop_iface : public objbase
    {
    public:
      virtual ~prop_iface ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:
      //! add new floating point property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_f (const fp_type_t def_value, const std::string  &short_name, const std::string  &description) = 0;
      
      //! add new integer property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_i (const i_type_t def_value, const std::string  &short_name, const std::string  &description) = 0; 

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_s (const s_type_t def_value, const std::string  &short_name, const std::string  &description) = 0; 

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_b (const b_type_t def_value, const std::string  &short_name, const std::string  &description) = 0; 
 
      //! clear all
      virtual void clear () = 0;

      //! return property value
      virtual const fp_type_t get_f (const std::string &name) const = 0;

      //! return property value
      virtual const i_type_t get_i (const std::string &name) const = 0;

      //! return property value
      virtual const s_type_t get_s (const std::string &name) const = 0;

      //! return property value
      virtual const b_type_t get_b (const std::string &name) const = 0;

      //! set value
      virtual void set_f (const std::string &name, const fp_type_t  value) = 0;

      //! set value
      virtual void set_i (const std::string &name, const i_type_t  value) = 0;

      //! set value
      virtual void set_s (const std::string &name, const s_type_t  value) = 0;

      //! set value
      virtual void set_b (const std::string &name, const b_type_t  value) = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_f (const std::string &name) const = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_i (const std::string &name) const = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_s (const std::string &name) const = 0;

      //! check (return false if property set by default, true otherwise
      virtual bool check_b (const std::string &name) const = 0;

      //! reset to default value
      virtual void reset_f (const std::string &name) = 0;

      //! reset to default value
      virtual void reset_i (const std::string &name) = 0;

      //! reset to default value
      virtual void reset_s (const std::string &name) = 0;

      //! reset to default value
      virtual void reset_b (const std::string &name) = 0;

      //! reset all
      virtual void reset_all () = 0;
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

    };
};

#endif //__PROP_IFACE_H


