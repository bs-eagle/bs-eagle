/** 
 * @file pvt_dead_oil.cpp
 * @brief realization of Dead Oil PVT properties
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-25
 */
#include "pvt_dead_oil.h"
using namespace std;
using namespace boost::python;


const std::string s_dens_idx = "surface_density";
const std::string m_dens_idx = "molar_density";
const std::string max_p_idx = "max_pressure";
const std::string min_p_idx = "min_pressure";
const std::string n_intervals_idx = "n_intervals";

namespace blue_sky
{
  pvt_dead_oil::pvt_dead_oil (bs_type_ctor_param) 
        : sp_in_prop (BS_KERNEL.create_object ("prop")),
        sp_in_table (BS_KERNEL.create_object ("table")),
        sp_calc_table (BS_KERNEL.create_object ("table"))

    {
      surface_density = molar_density = -1;
      init_prop ();
    }
  pvt_dead_oil::pvt_dead_oil (const pvt_dead_oil & /*src*/)
        : bs_refcounter (),
        sp_in_prop (BS_KERNEL.create_object ("prop")),
        sp_in_table (BS_KERNEL.create_object ("table")),
        sp_calc_table (BS_KERNEL.create_object ("table"))
    {
      surface_density = molar_density = -1;
      init_prop ();
    }

   void
   pvt_dead_oil::init_prop ()
    {
      sp_in_prop->add_property_f (-1.0, s_dens_idx, 
                                  std::string ("Surface density for phase"));
      sp_in_prop->add_property_f (-1.0, m_dens_idx, 
                                  std::string ("Molar density for phase"));
      sp_in_prop->add_property_f (1000.0, max_p_idx, 
                                  std::string ("Maximum allowed pressure for PVT table"));
      sp_in_prop->add_property_f (1.0, min_p_idx, 
                                  std::string ("Minimum allowed pressure for PVT table"));
      sp_in_prop->add_property_i (30, n_intervals_idx, 
                                  std::string ("Length of the calculated PVT table"));
    }

  int 
  pvt_dead_oil::check () const
    {
      // TODO: add input table checking
      return 0;
    }

  int 
  pvt_dead_oil::update ()
    {
      // TODO: add input table checking
      return 0;
    }


/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (pvt_dead_oil);
  BLUE_SKY_TYPE_STD_COPY (pvt_dead_oil);

  BLUE_SKY_TYPE_IMPL (pvt_dead_oil, pvt_oil_iface, "pvt_dead_oil", "Dead oil PVT properties", "Realization of PVT properties for Dead Oil");
}  // blue_sky namespace
