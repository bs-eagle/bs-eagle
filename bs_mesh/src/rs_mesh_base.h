#ifndef RS_MESH_BASE_H
#define RS_MESH_BASE_H
/*!
	\file rs_mesh_base.h
  \brief This file declares base class for reservoir simulation meshes
  \author Mark Khait
  \date 2009-07-20
 */

#include "flux_connections_iface.h"
#include "mesh_base.h"
#include "bs_hdf5_storage_v2.h"

using namespace blue_sky;

class BS_API_PLUGIN rs_mesh_base : public mesh_base
  {

  //-----------------------------------------
  // TYPES
  //-----------------------------------------

  public:
    ///////////////////////
    // BASE TYPES
    ///////////////////////
    typedef mesh_base base_t;

    typedef base_t::index_t                    index_t;

    typedef base_t::index_array_t              index_array_t;
    typedef base_t::item_array_t               item_array_t;

    typedef base_t::sp_bcsr_t                  sp_bcsr_t;
    typedef base_t::sp_idata_t                 sp_idata_t;
    
    ///////////////////////
    // OWN TYPES
    ///////////////////////

    typedef typename strategy_t::item_t                 item_t;
    
    typedef flux_connections_iface          flux_conn_iface_t;
    typedef smart_ptr <flux_conn_iface_t, true>         sp_flux_conn_iface_t;
    
    
  //-----------------------------------------
  //  METHODS
  //-----------------------------------------

  public:

    ///////////////////////
    // INIT
    ///////////////////////
        
    //! default constructor
    rs_mesh_base ()	{};

    //! default destructor
    virtual ~rs_mesh_base ()	{};
    
    //! init mesh
    void init_props (const sp_idata_t &idata);
    
    //! initialize int_to_ext indexation
    int init_int_to_ext();
    
    //! check mesh data
    void check_data () const;

    ///////////////////////
    // ACCESS
    ///////////////////////
    
    //! set darcy constant for correct transmissibility calculation 
    void set_darcy (double darcy_constant_)
    {
      darcy_constant = darcy_constant_;
    }
    
    //! return depths of cell centers (length n_active_elements)
    const item_array_t & get_depths () const
    {
      return depths;
    }
    
    ///////////////////////
    // WORK
    ///////////////////////

    //! allocate jacobian 
    virtual int build_jacobian_and_flux_connections (const sp_bcsr_t &/*jacobian*/, const sp_flux_conn_iface_t &/*flux_conn*/, 
                                                     index_array_t &/*boundary_array*/) = 0;

    // storage interface
    hdf5_group_v2 &
    save_info (hdf5_group_v2 &group) const;

    hdf5_group_v2 &
    save_data (hdf5_group_v2 &group) const;

  //-----------------------------------------
  //  VARIABLES
  //-----------------------------------------

  protected:
    
    item_t minpv;	  //!< minimum pore volume
    item_t minsv;	  //!< minimum volume for splicing
    item_t max_thickness;	  //!< maximum thickness between blocks for splicing
    
    auto_value <double> darcy_constant; //Darcy coefficient

    //smart_ptr on properties
    array_uint8_t   sp_actnum;	//!< smart_ptr on actnum array
    array_float16_t sp_poro;		//!< smart_ptr on poro array
    array_float16_t sp_ntg; 		//!< smart_ptr on ntg array
    array_float16_t sp_multpv;	//!< smart_ptr on multpv array

    item_array_t depths; //!< depths of elements center points

  };
#endif // RS_MESH_BASE_H
