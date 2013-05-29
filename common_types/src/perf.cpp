/** 
 * @file perf.cpp
 * @brief implementation of perf storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include <iomanip>
#include <iostream>


#include "bs_kernel.h"
#include "perf.h"


using namespace boost;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{

  perf::perf (bs_type_ctor_param) 
    {
      sp_prop = BS_KERNEL.create_object ("prop");
      if (!sp_prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      
    }
  perf::perf (const perf& rhs) 
        : bs_refcounter ()
    {
      *this = rhs;
    }
  void
  perf::init_prop ()
    {
      sp_prop->add_property_f (-1, L"parent_md", 
                               L"Position of the fracture center in the well branch (m.)");
      sp_prop->add_property_s (L"main", L"parent", 
                               L"Parent well branch (m.)");
      sp_prop->add_property_f (1, L"length", 
                               L"Length (m.) of perforation interval");
      sp_prop->add_property_f (0.06, L"radius", 
                               L"Radius (m.) of the wellbore");
      sp_prop->add_property_f (-1, L"trans", 
                               L"Transmissibility (cP-m3/day-bars) of the perforation");
      sp_prop->add_property_f (-1, L"kh", 
                               L"Effective Kh (mD-m) of the perforation");
      sp_prop->add_property_f (0, L"skin", 
                               L"Skin factor of the perforation");
      sp_prop->add_property_f (-1, L"r0", 
                               L"Pressure equivalent radius");

    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  perf::py_str () const
    {
      std::stringstream s;
      s << "\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (perf);
  BLUE_SKY_TYPE_STD_COPY (perf);

  BLUE_SKY_TYPE_IMPL(perf,  perf_iface, "perf", "perf storage", "realization of well perf storage");

}  // blue_sky namespace
