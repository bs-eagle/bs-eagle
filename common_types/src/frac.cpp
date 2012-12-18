/** 
 * @file frac.cpp
 * @brief implementation of frac storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include <iomanip>
#include <iostream>


#include "bs_kernel.h"
#include "frac.h"


using namespace boost;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{

  frac::frac (bs_type_ctor_param) 
    {
      sp_prop = BS_KERNEL.create_object ("prop");
      if (!sp_prop)
        {
          bs_throw_exception ("Type (prop) not refractered");
        }
      
    }
  frac::frac (const frac& rhs) 
        : bs_refcounter ()
    {
      *this = rhs;
    }
  void
  frac::init_prop ()
    {
      sp_prop->add_property_b (true, "is_vertical", 
                               L"Should be true if fracture is vertical");
      sp_prop->add_property_b (true, "is_symmetric", 
                               L"Is fracture is symmetrical");
      sp_prop->add_property_b (true, "inf_perm", 
                               L"If true use infinum permeability for fracture");
      sp_prop->add_property_f (100000.0, "perm", 
                               L"Permiability (mD) for fracture (only if inf_perm == false)");
      sp_prop->add_property_f (0.005, "wf", 
                               L"Width of the fracture (m.) (only if inf_perm == false)");
      sp_prop->add_property_f (50, "half_length", 
                               L"Half length of the fracture (m.)");
      sp_prop->add_property_f (20, "up_half_height", 
                               L"Half height in up direction of the fracture (m.)");
      sp_prop->add_property_f (20, "down_half_height", 
                               L"Half height in down direction of the fracture (m.)");
      sp_prop->add_property_f (0, "angle", 
                               L"Fracture angle (grad)");
      sp_prop->add_property_f (50, "hor_main_radius", 
                               L"Main radius for horizontal well (m.)");
      sp_prop->add_property_f (50, "hor_sec_radius", 
                               L"Secondary radius for horizontal well (m.)");
      sp_prop->add_property_f (-1, "parent_md", 
                               L"Position of the fracture center in the well branch (m.)");
      sp_prop->add_property_s ("main", "parent", 
                               L"Parent well branch (m.)");

    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  frac::py_str () const
    {
      std::stringstream s;
      s << sp_prop->py_str () << "\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (frac);
  BLUE_SKY_TYPE_STD_COPY (frac);

  BLUE_SKY_TYPE_IMPL(frac,  frac_iface, "frac", "frac storage", "realization of well frac storage");

}  // blue_sky namespace
