#ifndef BS_MESH_GRDECL_H
#define BS_MESH_GRDECL_H
/*!
	\file bs_mesh_grdecl.h
  \brief This file declares bs wrapper over mesh_grdecl
  \author Mark Khait
	\date 2008-05-20
 */

#include "mesh_grdecl.h"
#include "rs_smesh_iface.h"


namespace blue_sky
  {

  class idata;

  template <typename strategy_t>
  class BS_API_PLUGIN mesh_grdecl_keywords;

  template<class strategy_t>
  class BS_API_PLUGIN bs_mesh_grdecl : virtual public rs_smesh_iface<strategy_t>
    {

    ///////////////////////////////
    //  INTERNAL TYPE DECLARATION
    ///////////////////////////////
    public:
      ///////////////////////
      // BASE TYPES
      ///////////////////////
      typedef rs_smesh_iface<strategy_t>                  base_t;

      typedef typename base_t::index_t                    index_t;
      typedef typename base_t::item_t                     item_t;

      typedef typename base_t::index_array_t              index_array_t;
      typedef typename base_t::item_array_t               item_array_t;

      typedef typename base_t::sp_flux_conn_iface_t       sp_flux_conn_iface_t;
      typedef typename base_t::sp_bcsr_t                  sp_bcsr_t;
      typedef typename base_t::sp_idata_t                 sp_idata_t;
      typedef typename base_t::point3d_t                  point3d_t;

      ///////////////////////
      // OWN TYPES
      ///////////////////////

    public:
      //! blue-sky class declaration
      BLUE_SKY_TYPE_DECL(bs_mesh_grdecl);

      //! default destructor
      ~bs_mesh_grdecl ()	{};

      ///////////////////////
      // INITIALIZATION
      ///////////////////////

      //! init mesh properties
      void init_props (const sp_idata_t &idata)
        {wrapped.init_props (idata);};

      //! initialize int_to_ext indexation
      int init_int_to_ext()
        {return wrapped.init_int_to_ext ();};

      //! initialize int_to_ext indexation
      int init_ext_to_int()
        {return wrapped.init_ext_to_int();};


      ///////////////////////
      // ACCESS VARIABLES
      ///////////////////////


      //! return number of active mesh elements
      index_t get_n_active_elements () const
        {return wrapped.get_n_active_elements ();};

      //! return number of mesh elements
      index_t get_n_elements () const
        {return wrapped.get_n_elements ();};

      //! return number of mesh elements connections
      index_t get_n_connections () const
        {return wrapped.get_n_connections ();};

      //! return mesh dimensions range
      void get_dimensions_range (item_t &dim1_max, item_t &dim1_min,
        item_t &dim2_max, item_t &dim2_min,
        item_t &dim3_max, item_t &dim3_min) const
        {return wrapped.get_min_max_xyz (dim1_max, dim1_min, dim2_max, dim2_min, dim3_max, dim3_min);};

      //! get mesh dimensions
      typename base_t::index_point3d_t get_dimens ()
      {return wrapped.get_dimens();};

      //! return element size in all 3 dimensions
      void get_element_size (const index_t n_element, item_t &d_dim1, item_t &d_dim2, item_t &d_dim3) const
        {wrapped.get_block_dx_dy_dz(n_element, d_dim1, d_dim2, d_dim3);};

      //! return element size in 3rd dimension
      item_t get_element_dim3_size (const index_t n_element) const
        {return wrapped.get_block_dz(n_element);};

      //! return center point of an element
      point3d_t get_element_center (const index_t i, const index_t j, const index_t k)const
        {return wrapped.get_center(i, j, k);};

      //! return center point of an element by internal element number
      point3d_t get_element_center (const index_t n_element)const
        {return wrapped.get_center(n_element);};

      //! return depth of mesh element
      item_t get_element_depth(const index_t n_element) const
        {return wrapped.get_depth(n_element);};

      //! return number of mesh elements connections
      item_t get_element_dtop(const index_t n_element) const
        {return wrapped.get_dtop(n_element);};

      //! get element internal number by external
      index_t get_element_ext_to_int (const index_t n_element) const
        {return wrapped.get_element_ext_to_int(n_element);};

      //! get element external number by internal
      index_t get_element_int_to_ext (const index_t n_element) const
        {return wrapped.get_element_int_to_ext(n_element);};

      //! return I, J and K structured mesh coordinates of an element by internal number
      void get_element_int_to_ijk (const index_t n_element, index_t &i, index_t &j, index_t &k) const
        {wrapped.inside_to_XYZ(n_element, i, j, k);};

      //! return internal number of an element by I, J and K structured mesh coordinates
      index_t get_element_ijk_to_int (const index_t i, const index_t j, const index_t k) const
        {return wrapped.XYZ_to_inside(i, j, k);};

      //! return coords of block vertexes by IJK indexes
      grd_ecl::fpoint3d_vector top_cube (const index_t i, const index_t j, const index_t k) const
        {return wrapped.top_cube (i, j, k);};

      //! return coords of block vertexes by n_block index
      grd_ecl::fpoint3d_vector top_cube (const index_t index) const
        {return wrapped.top_cube (index);};

      ///////////////////////
      // ACCESS ARRAYS
      ///////////////////////


      //! get const int_to_ext
      const index_array_t & get_int_to_ext() const
        {return wrapped.get_int_to_ext ();};

      //! get const ext_to_int
      const index_array_t & get_ext_to_int() const
        {return wrapped.get_ext_to_int ();};

      //! get mesh elements volumes
      const item_array_t &get_volumes () const
        {return wrapped.get_volumes ();};


      //! return depths of cell centers (length n_active_elements)
      const item_array_t &get_depths () const
        {return wrapped.get_depths();};

      //! set darcy constant for correct transmissibility calculation
      void set_darcy (double darcy_constant_)
        {wrapped.set_darcy (darcy_constant_);};

      ///////////////////////
      // WORK
      ///////////////////////

      //! check data
      void check_data () const
        {wrapped.check_data();};

      //! allocate jacobian
      int build_jacobian_and_flux_connections (const sp_bcsr_t &jacobian, const sp_flux_conn_iface_t &flux_conn, index_array_t &boundary_array)
        {return wrapped.build_jacobian_and_flux_connections (jacobian, flux_conn, boundary_array);};


    ////////////////////
    // wrapped class
    ///////////////////

    private:

      mesh_grdecl<strategy_t> wrapped;
    };

};//namespace blue_sky
#endif // BS_MESH_GRDECL_H
