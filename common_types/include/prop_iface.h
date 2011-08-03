/** 
 * @file prop_iface.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-27
 */
#ifndef __PROP_IFACE_H
#define __PROP_IFACE_H
#include <string>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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
      typedef boost::archive::text_iarchive           tia_t;
      typedef boost::archive::text_oarchive           toa_t;

      virtual ~prop_iface ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:
      //! add new floating point property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_f (const t_double def_value, const std::string  &short_name, const std::string  &description) = 0;
      
      //! add new integer property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_i (const t_long def_value, const std::string  &short_name, const std::string  &description) = 0; 

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_s (const std::string def_value, const std::string  &short_name, const std::string  &description) = 0; 

      //! add new string property to the list (return index of the property or < 0 if error occur)
      virtual void add_property_b (const bool def_value, const std::string  &short_name, const std::string  &description) = 0; 
 
      //! clear all
      virtual void clear () = 0;

      //! return property value
      virtual t_double get_f (const std::string &name) const = 0;

      //! return property value
      virtual t_long get_i (const std::string &name) const = 0;

      //! return property value
      virtual std::string get_s (const std::string &name) const = 0;

      //! return property value
      virtual bool get_b (const std::string &name) const = 0;

      //! set value
      virtual void set_f (const std::string &name, const t_double  value) = 0;

      //! set value
      virtual void set_i (const std::string &name, const t_long  value) = 0;

      //! set value
      virtual void set_s (const std::string &name, const std::string  value) = 0;

      //! set value
      virtual void set_b (const std::string &name, const bool  value) = 0;

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

      virtual void save (toa_t &ar) const = 0;
      virtual void load (tia_t &ar) = 0;
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

    };
};

#endif //__PROP_IFACE_H


