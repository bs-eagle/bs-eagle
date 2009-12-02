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

  idata::~idata ()
  {

  }

  idata::idata(bs_type_ctor_param param)
  : bs_node(bs_node::create_node (new this_t::idata_traits ())),
  rpo_model (0), // RPO_DEFAULT_MODEL
  minimal_pore_volume (DEFAULT_MINIMAL_PORE_VOLUME),
  minimal_splice_volume (DEFAULT_MINIMAL_SPLICE_VOLUME),
  maximum_splice_thickness (DEFAULT_MAXIMUM_SPLICE_THICKNESS),
  pvt_region (1),
  sat_region (1),
  eql_region (1),
  fip_region (1),
  fi_n_phase (0),
  fi_phases (0),
  rock_region (1),
  i_map (BS_KERNEL.create_object (amap_i::bs_type())),
  d_map (BS_KERNEL.create_object (amap_d::bs_type())),
  init_section (0)
  {
    init();
  }

  idata::idata(const this_t &src)
      : bs_refcounter (src), bs_node(src), 
      i_map(give_kernel::Instance().create_object_copy(src.i_map)),
      d_map(give_kernel::Instance().create_object_copy(src.d_map))
      //scal3(give_kernel::Instance().create_object_copy(src.scal3)),
  {
    *this = src;
  }

  void idata::init()
  {
    //depth.resize((nx+1) * (ny+1) * (nz+1));
		ahelper.init_names_maps ();
  }

  idata &idata::operator=(const this_t &src)
  {
    i_map = src.i_map;
    d_map = src.d_map;
		ahelper = src.ahelper;
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

  void idata::set_defaults_in_pool()
  {
    //d_map->create_item (MULTX,  &d_pool_sizes[ARRAY_POOL_TOTAL * MULTX],  d_pool_default_values[MULTX]);
    //d_map->create_item (MULTY,  &d_pool_sizes[ARRAY_POOL_TOTAL * MULTY],  d_pool_default_values[MULTY]);
    //d_map->create_item (MULTZ,  &d_pool_sizes[ARRAY_POOL_TOTAL * MULTZ],  d_pool_default_values[MULTZ]);
    //d_map->create_item (NTG,    &d_pool_sizes[ARRAY_POOL_TOTAL * NTG],    d_pool_default_values[NTG]);
    //d_map->create_item (MULTPV, &d_pool_sizes[ARRAY_POOL_TOTAL * MULTPV], d_pool_default_values[MULTPV]);
    i_map->create_item (ACTNUM, &i_pool_sizes[ARRAY_POOL_TOTAL * ACTNUM], i_pool_default_values[ACTNUM]);
  }


  void idata::set_region (int r_pvt,int r_sat, int r_eql, int r_fip)
  {
    this->pvt_region = r_pvt;
    this->fip_region = r_fip;
    this->sat_region = r_sat;
    this->eql_region = r_eql;

    // check
    if (r_pvt <= 0 || r_sat <= 0 || r_eql <= 0 || r_fip <= 0)
      {
        bs_throw_exception ("One of init parameters <= 0");
      }

    assign (this->rock, r_pvt, -1);
    assign (this->p_ref, r_pvt, -1);

    this->equil.resize(EQUIL_TOTAL * eql_region); //!TODO: EQUIL_TOTAL instead of 3

    this->pvto.resize(r_pvt);
    this->pvtdo.resize(r_pvt);
    this->pvtg.resize(r_pvt);
    this->pvtw.resize(r_pvt);
  }

  void idata::set_density (const shared_vector <double> &density)
  {
    if ((density.size() % 3 != 0) || (density.size()<3))
      {
        bs_throw_exception ("Not enough valid arguments yet");
      }

    set_density_internal(&density[0]);
  }

  void idata::set_density_internal (const double *density)
  {
    std::ostringstream out_s;

    if (!pvto.size())
      throw bs_exception("idata.set_density()","PVT table for oil has not been initialized yet");

    if (!pvtw.size())
      throw bs_exception("idata.set_density()","PVT table for water has not been initialized yet");

    if (!pvtg.size())
      throw bs_exception("idata.set_density()","PVT table for gas has not been initialized yet");

    if (!rock.size())
      throw bs_exception("idata.set_density()","Rock properties table has not been initialized yet");

    if (!equil.size())
      throw bs_exception("idata.set_density()","EQUIL table has not been initialized yet");

    for (size_t i = 0; i < pvt_region; ++i)
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
  
	/* 
			ARRAYS_HELPER methods
	*/

	idata::arrays_helper::arrays_helper () {
		dummy_array_i = new uint8_t [1];
		dummy_array_d = new float16_t [1];
	}

	idata::arrays_helper::~arrays_helper () {
		delete [] dummy_array_i;
		delete [] dummy_array_d;
	}

	void idata::arrays_helper::init_names_maps () {
		for (int i = MPFANUM; i < ARR_TOTAL_INT; ++i)
			add_correspondence_i(int_names_table[i], i);

		for (int i = SGL; i < ARR_TOTAL_DOUBLE; ++i)
		add_correspondence_d(double_names_table[i], i);
	}

	void idata::arrays_helper::add_correspondence_i (const std::string &name, int index) {
		names_map_i[name] = index;
	}
	void idata::arrays_helper::add_correspondence_d (const std::string &name, int index) {
		names_map_d[name] = index;
	}

	int idata::arrays_helper::get_idx_i (const std::string &name) const {
		names_map_t::const_iterator iter = names_map_i.find (name);
		if (iter != names_map_i.end ())
			return iter->second;
		return -1;
	}
	int idata::arrays_helper::get_idx_d (const std::string &name) const {
		names_map_t::const_iterator iter = names_map_d.find (name);
		if (iter != names_map_d.end ())
			return iter->second;
		return -1;
	}

  array_uint8_t
  idata::get_int_non_empty_array (int array_index)
  {
    return tools::get_non_empty ((*i_map)[array_index].array);
  }
  array_float16_t
  idata::get_float_non_empty_array (int array_index)
  {
    return tools::get_non_empty ((*d_map)[array_index].array);
  }
  
  array_uint8_t
  idata::get_int_array (const std::string & array_name)
  {
    int array_index = ahelper.get_idx_i(array_name);
    array_uint8_t dummy;
    if (array_index < 0)
      {
        return dummy;
      }
    else
      {
        return (*i_map)[array_index].array;
      }
  }
  array_float16_t
  idata::get_float_array (const std::string &array_name)
  {
    int array_index = ahelper.get_idx_d(array_name);
    array_float16_t dummy;
    if (array_index < 0)
      {
        return dummy;
      }
    else
      {
        return (*d_map)[array_index].array;
      }
  }
  
  array_uint8_t &
  idata::get_int_non_empty_array (const std::string &array_name)
  {
    int array_index = ahelper.get_idx_i(array_name);
    if (array_index < 0)
      {
        bs_throw_exception (boost::format ("Integer array %s is not initialized yet!") % array_name);
      }
    else
      {
        return (*i_map)[array_index].array;
      }
  }
  
  array_float16_t
  idata::get_float_non_empty_array (const std::string &array_name)
  {
    int array_index = ahelper.get_idx_d(array_name);
    if (array_index < 0)
      {
        bs_throw_exception (boost::format ("Float array %s is not initialized yet!") % array_name);
      }
    else
      {
        return (*d_map)[array_index].array;
      }
  }
  
  bool
  idata::contain (const std::string &array_name) const
  {
    int array_index = ahelper.get_idx_i(array_name);
    if (array_index >= 0) return true;
    
    array_index = ahelper.get_idx_d(array_name);
    if (array_index >= 0) return true;
    
    return false;    
  }

  // create object
  BLUE_SKY_TYPE_STD_CREATE (idata);
  BLUE_SKY_TYPE_STD_COPY (idata);

  BLUE_SKY_TYPE_IMPL (idata, bs_node, "idata", "Initial data storage", "Initial data storage");

  /*
  	\brief register plugin
   */
  bool register_idata (const plugin_descriptor& pd)
  {
    bool res;
    res = BLUE_SKY_REGISTER_TYPE (pd, idata); 
    return res;
  }
}
