/**
	\file data_class.cpp
	\brief implimenataion of idata class methods
	\author Nikonov Maxim
*/
#include "bs_bos_core_data_storage_stdafx.h"

#include "data_class.h"
#include "arrays.h"
#include "arrays_tables.h"
#include "constants.h"

//! Default minimal pore volume allowed for active cells
//#define DEFAULT_MINIMAL_PORE_VOLUME     1.e-1
#define DEFAULT_MINIMAL_PORE_VOLUME     1.e-7
//! Default minimal pore volume allowed for active cells to splice with other cells
// miryanov: set MINSV equal to MINPV for SPE10 model
#define DEFAULT_MINIMAL_SPLICE_VOLUME     1.e-7
//! Default maximum thickness allowed between active cells to be coupled
#define DEFAULT_MAXIMUM_SPLICE_THICKNESS     1.e-2


namespace blue_sky
  {
  /*!
  	\class idata_traits
  	\brief for node insertion control
  */
  
  struct idata::idata_traits : public bs_node::sort_traits
    {
      struct idata_key : bs_node::sort_traits::key_type
        {
          virtual bool sort_order(const key_ptr&) const
            {
              return true;
            }
        };

      virtual const char *sort_name() const
        {
          return "idata trait";
        }

      virtual key_ptr key_generator(const sp_link&) const
        {
          return new idata_key();
        }

      virtual bool accepts(const sp_link &l) const
        {
          if (l->name().find("_in_pool",0,l->name().size()))
            return smart_ptr<this_t>(l->data(), bs_dynamic_cast());
          return false;
        }
    };

  idata::pvt_info::pvt_info ()
  : main_data_ (BS_KERNEL.create_object (v_float::bs_type ()))
  {
  }

  idata::~idata ()
  {
    h5_pool->close_file();
  }

  
  idata::idata(bs_type_ctor_param /*param*/)
  : bs_node(bs_node::create_node (new this_t::idata_traits ())),
  props (BS_KERNEL.create_object ("prop")),
  h5_pool (BS_KERNEL.create_object ("h5_pool")),
  equil_regions (BS_KERNEL.create_object (v_int::bs_type ())),
  rock (BS_KERNEL.create_object (v_float::bs_type ())),
  p_ref (BS_KERNEL.create_object (v_float::bs_type ())),
  equil (BS_KERNEL.create_object (v_float::bs_type ()))
  {
    init();
  }

  
  idata::idata(const this_t &src)
  : bs_refcounter (src), bs_node(src), 
  props(give_kernel::Instance().create_object_copy(src.props)),
  h5_pool(give_kernel::Instance().create_object_copy(src.h5_pool))
  //equil_regions (BS_KERNEL.create_object_copy (src.equil_regions)),
  //rock (BS_KERNEL.create_object_copy (src.rock)),
  //p_ref (BS_KERNEL.create_object_copy (src.p_ref)),
  //equil (BS_KERNEL.create_object_copy (src.equil))
  {
    *this = src;
  }

  
  void idata::init()
  {
    //depth.resize((nx+1) * (ny+1) * (nz+1));
    h5_pool->open_file ("bs_data_storage.h5", "/pool");
    
    props->add_property_i (0, "rpo_model", "3-ph oil relative permeability model: flag 0, 1 or 2 (stone model)");  
    props->add_property_i (1, "pvt_region", "Number of PVT regions in simulation");
    props->add_property_i (1, "sat_region", "Number of saturation regions in simulation");
    props->add_property_i (1, "eql_region", "Number of equilibrium regions in simulation");
    props->add_property_i (1, "fip_region", "Number of FIP regions in simulation");
    props->add_property_i (1, "rock_region", "Number of ROCK regions");
    props->add_property_i (0, "init_section", "flag indicating whether we have init section");
    
    props->add_property_s ("", "title", "Current model title");
    props->add_property_b (0, "oil_phase", "True if oil phase exists");
    props->add_property_b (0, "water_phase", "True if water phase exists");
    props->add_property_b (0, "gas_phase", "True if gas phase exists");
    props->add_property_b (0, "scalecrs", "True if SCALECRS is enabled");
  }

  void idata::flush_pool()
  {
    h5_pool->flush ();
  }

  
  idata &idata::operator=(const this_t &src)
  {
    h5_pool = src.h5_pool;
    return *this;
  }


  
  idata::vval_vs_depth &idata::get_prvd()
  {
    return prvd;
  }

  
  idata::vval_vs_depth &idata::get_rsvd()
  {
    return rsvd;
  }

  
  idata::vval_vs_depth &idata::get_pbvd()
  {
    return pbvd;
  }

  //idata::vsp_pvt &idata::get_pvt()
  //{
  //  return pvt;
  //}

  //const idata::sp_scal3p &idata::get_scal()
  //{
  //  return scal3;
  //}
  /*!
  	\brief updating dx, dy, dz
  				 add values in nonactive blocks for mesh generation algorithm
  	\return if success  0
  */

  
  void idata::set_region (int r_pvt,int r_sat, int r_eql, int r_fip)
  {
    props->set_i ("pvt_region", r_pvt);
    props->set_i ("fip_region", r_fip);
    props->set_i ("sat_region", r_sat);
    props->set_i ("eql_region", r_eql);

    // check
    if (r_pvt <= 0 || r_sat <= 0 || r_eql <= 0 || r_fip <= 0)
      {
        bs_throw_exception ("One of init parameters <= 0");
      }
    
    t_long def_val = -1;
    
    rock->init (r_pvt, def_val);
    p_ref->init (r_pvt, def_val);

    equil->resize (EQUIL_TOTAL * r_eql); //!TODO: EQUIL_TOTAL instead of 3

    pvto.resize (r_pvt);
    pvtdo.resize (r_pvt);
    pvtg.resize (r_pvt);
    pvtw.resize (r_pvt);
  }

  
  void idata::set_density (spv_float density)
  {
    if ((density->size() % 3 != 0) || (density->size()<3))
      {
        bs_throw_exception ("Not enough valid arguments yet");
      }

    set_density_internal(&(*density)[0]);
  }

  
  void idata::set_density_internal (const t_float *density)
  {
    std::ostringstream out_s;

    if (!pvto.size())
      throw bs_exception("idata.set_density()","PVT table for oil has not been initialized yet");

    if (!pvtw.size())
      throw bs_exception("idata.set_density()","PVT table for water has not been initialized yet");

    if (!pvtg.size())
      throw bs_exception("idata.set_density()","PVT table for gas has not been initialized yet");

    if (!rock->size())
      throw bs_exception("idata.set_density()","Rock properties table has not been initialized yet");

    if (!equil->size())
      throw bs_exception("idata.set_density()","EQUIL table has not been initialized yet");

    t_int n_pvt = props->get_i ("pvt_region");

    for (t_int i = 0; i < n_pvt; ++i)
      {
        idata::pvt_info &pvto__ = pvto[i];
        idata::pvt_info &pvtw__ = pvtw[i];
        idata::pvt_info &pvtg__ = pvtg[i];
        idata::pvt_info &pvdo__ = pvtdo[i];

        pvto__.has_density_		= true;
        pvto__.density_				= density[i*3];
        pvto__.molar_density_ = density[i*3];

        pvdo__.has_density_		= true;
        pvdo__.density_				= density[i*3];
        pvdo__.molar_density_ = density[i*3];

        pvtw__.has_density_		= true;
        pvtw__.density_				= density[i*3+1];
        pvtw__.molar_density_ = density[i*3+1];

        pvtg__.has_density_		= true;
        pvtg__.density_				= density[i*3+2];
        pvtg__.molar_density_ = density[i*3+2];
      }
  }
  
  
  spv_int idata::get_i_array (const std::string & array_name, bool safe)
  {
    return safe 
        ? h5_pool->get_i_data (array_name)
        : h5_pool->get_i_data_unsafe (array_name)
        ;
  }
  
  spv_float idata::get_fp_array (const std::string &array_name, bool safe)
  {
    return safe
        ? h5_pool->get_fp_data (array_name)
        : h5_pool->get_fp_data_unsafe (array_name)
        ;
  }

  bool
  idata::contains_i_array (const std::string &array_name) 
  {
    const spv_int &a = h5_pool->get_i_data_unsafe (array_name);
    return a && a->size ();
  }

  bool
  idata::contains_fp_array (const std::string &array_name)
  {
    const spv_float &a = h5_pool->get_fp_data_unsafe (array_name);
    return a && a->size ();
  }
  
  int idata::set_i_array (const std::string & array_name,  spv_int array)
  {
    return h5_pool->set_i_data (array_name, array);
  }
  

  int idata::set_fp_array (const std::string & array_name,  spv_float array)
  {
    return h5_pool->set_fp_data (array_name, array);
  }

  spv_int idata::create_i_array (const std::string & /*array_name*/,  t_int *array_dimens, t_int def_value)
  {
    spv_int new_array;
    t_long n;
    npy_intp dims[3];

    dims[0] = array_dimens[ARRAY_POOL_NX_A] * props->get_i ("nx") + array_dimens[ARRAY_POOL_NX_B];
    dims[1] = array_dimens[ARRAY_POOL_NY_A] * props->get_i ("ny") + array_dimens[ARRAY_POOL_NY_B];
    dims[2] = array_dimens[ARRAY_POOL_NZ_A] * props->get_i ("nz") + array_dimens[ARRAY_POOL_NZ_B];
    n = dims[0] * dims[1] * dims[2];

    new_array = BS_KERNEL.create_object (v_int::bs_type ());
    new_array->resize (n);
    std::fill (new_array->begin (), new_array->end(), def_value);
    new_array->reshape (3, dims);

    return new_array;
  }
  

  spv_float idata::create_fp_array (const std::string & /*array_name*/,  t_int *array_dimens, t_float def_value)
  {
    spv_float new_array;
    t_long n;
    npy_intp dims[3];
    dims[0] = array_dimens[ARRAY_POOL_NX_A] * props->get_i ("nx") + array_dimens[ARRAY_POOL_NX_B];
    dims[1] = array_dimens[ARRAY_POOL_NY_A] * props->get_i ("ny") + array_dimens[ARRAY_POOL_NY_B];
    dims[2] = array_dimens[ARRAY_POOL_NZ_A] * props->get_i ("nz") + array_dimens[ARRAY_POOL_NZ_B];
    n = dims[0] * dims[1] * dims[2];

    new_array = BS_KERNEL.create_object (v_float::bs_type ());
    new_array->resize (n);
    std::fill (new_array->begin (), new_array->end(), def_value);
    new_array->reshape (3, dims);

    return new_array;
  }

  // create object
  BLUE_SKY_TYPE_STD_CREATE(idata)
  BLUE_SKY_TYPE_STD_COPY(idata)

  BLUE_SKY_TYPE_IMPL(idata, bs_node, "idata", "BOS_Core Initial data storage", "BOS_Core Initial data storage")
}
