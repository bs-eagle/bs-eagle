#ifndef RS_MESH_BASE_H
#define RS_MESH_BASE_H
/*!
	\file rs_mesh_base.h
  \brief This file declares base class for reservoir simulation meshes
  \author Mark Khait
  \date 2009-07-20
 */


#ifndef PURE_MESH
  #include "flux_connections_iface.h"
  #include "well_pool_iface.h"

  using namespace blue_sky;
#else
  #include "pure_mesh.h"
#endif

  #include "mesh_base.h"


class BS_API_PLUGIN rs_mesh_base : public mesh_base
  {

  //-----------------------------------------
  // TYPES
  //-----------------------------------------

  public:
    ///////////////////////
    // BASE TYPES
    ///////////////////////
    typedef mesh_base                       base_t;

    ///////////////////////
    // OWN TYPES
    ///////////////////////

#ifndef PURE_MESH    
    typedef flux_connections_iface                      flux_conn_iface_t;
    typedef smart_ptr <flux_conn_iface_t, true>         sp_flux_conn_iface_t;
    
    typedef BS_SP(well_pool_iface)                      sp_well_pool_t;
#else
    typedef csr_matrix                                  flux_conn_iface_t;
#endif
    
    
  //-----------------------------------------
  //  METHODS
  //-----------------------------------------

  public:

    ///////////////////////
    // INIT
    ///////////////////////
        
    //! default constructor
    rs_mesh_base ();

    //! default destructor
    virtual ~rs_mesh_base ()	{};

    //! init mesh from pool
    void init_props (const sp_hdm_t hdm);
    
    /*   
    //! init mesh from arrays
    void init_props (array_uint8_t   actnum_array,
                     array_float16_t poro_array,
                     array_float16_t ntg_array,
                     array_float16_t multpv_array);
    */                   
    
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
    const spv_float get_depths () const
    {
      return depths;
    }
    
    ///////////////////////
    // WORK
    ///////////////////////

    //! allocate jacobian 
    virtual int build_jacobian_and_flux_connections (const sp_bcsr_t /*jacobian*/, const sp_flux_conn_iface_t/*flux_conn*/, 
                                                     spv_long /*boundary_array*/) = 0;

    //! find well`s trajectories and mesh cells intersection
    virtual int intersect_trajectories () = 0;

  //-----------------------------------------
  //  VARIABLES
  //-----------------------------------------

  protected:

#ifndef PURE_MESH
    sp_well_pool_t well_pool; //!< sql well pool
#endif

    t_double minpv;	  //!< minimum pore volume
    t_double minsv;	  //!< minimum volume for splicing
    t_double max_thickness;	  //!< maximum thickness between blocks for splicing
    
    t_double darcy_constant; //Darcy coefficient

    spv_int   actnum_array;	//!< smart_ptr on actnum array
    spv_float poro_array;		//!< smart_ptr on poro array
    spv_float ntg_array; 		//!< smart_ptr on ntg array
    spv_float multpv_array;	//!< smart_ptr on multpv_array array
    spv_float depths;       //!< depths of elements center points

  };
#endif // RS_MESH_BASE_H
