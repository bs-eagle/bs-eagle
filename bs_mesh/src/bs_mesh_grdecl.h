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

  template <typename strategy_t>
  class idata;

  template <typename strategy_t>
  class mesh_grdecl_keywords;

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

      typedef typename base_t::i_type_t                   i_type_t;
      typedef typename base_t::fp_type_t                  fp_type_t;

      typedef typename base_t::sp_i_array_t               sp_i_array_t;
      typedef typename base_t::sp_fp_array_t              sp_fp_array_t;

      typedef typename base_t::sp_flux_conn_iface_t       sp_flux_conn_iface_t;
      typedef typename base_t::sp_bcsr_t                  sp_bcsr_t;
      typedef typename base_t::sp_idata_t                 sp_idata_t;
      typedef typename base_t::point3d_t                  point3d_t;

	  typedef bs_array< typename strategy_t::fp_storage_type_t > fp_storage_array_t;
	  typedef smart_ptr< fp_storage_array_t >                sp_fp_storage_array_t;

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
      i_type_t get_n_active_elements () const
        {return wrapped.get_n_active_elements ();};

      //! return number of mesh elements
      i_type_t get_n_elements () const
        {return wrapped.get_n_elements ();};

      //! return number of mesh elements connections
      i_type_t get_n_connections () const
        {return wrapped.get_n_connections ();};

      //! return mesh dimensions range
      void get_dimensions_range (fp_type_t &dim1_max, fp_type_t &dim1_min,
        fp_type_t &dim2_max, fp_type_t &dim2_min,
        fp_type_t &dim3_max, fp_type_t &dim3_min) const
        {return wrapped.get_min_max_xyz (dim1_max, dim1_min, dim2_max, dim2_min, dim3_max, dim3_min);};

      //! get mesh dimensions
      typename base_t::index_point3d_t get_dimens ()
      {return wrapped.get_dimens();};

      //! return element size in all 3 dimensions
      void get_element_size (const i_type_t n_element, fp_type_t &d_dim1, fp_type_t &d_dim2, fp_type_t &d_dim3) const
        {wrapped.get_block_dx_dy_dz(n_element, d_dim1, d_dim2, d_dim3);};

      //! return element size in 3rd dimension
      fp_type_t get_element_dim3_size (const i_type_t n_element) const
        {return wrapped.get_block_dz(n_element);};

      //! return center point of an element
      point3d_t get_element_center (const i_type_t i, const i_type_t j, const i_type_t k)const
        {return wrapped.get_center(i, j, k);};

      //! return center point of an element by internal element number
      point3d_t get_element_center (const i_type_t n_element)const
        {return wrapped.get_center(n_element);};

      //! return depth of mesh element
      fp_type_t get_element_depth(const i_type_t n_element) const
        {return wrapped.get_depth(n_element);};

      //! return number of mesh elements connections
      fp_type_t get_element_dtop(const i_type_t n_element) const
        {return wrapped.get_dtop(n_element);};

      //! get element internal number by external
      i_type_t convert_ext_to_int (const i_type_t n_element) const
        {return wrapped.convert_ext_to_int(n_element);};

      //! get element external number by internal
      i_type_t get_element_int_to_ext (const i_type_t n_element) const
        {return wrapped.get_element_int_to_ext(n_element);};

      //! return I, J and K structured mesh coordinates of an element by internal number
      void get_element_int_to_ijk (const i_type_t n_element, i_type_t &i, i_type_t &j, i_type_t &k) const
        {wrapped.inside_to_XYZ(n_element, i, j, k);};

      //! return internal number of an element by I, J and K structured mesh coordinates
      i_type_t get_element_ijk_to_int (const i_type_t i, const i_type_t j, const i_type_t k) const
        {return wrapped.XYZ_to_inside(i, j, k);};

      //! return coords of block vertexes by IJK indexes
      grd_ecl::fpoint3d_vector calc_element (const i_type_t i, const i_type_t j, const i_type_t k) const
        { typename mesh_grdecl<strategy_t>::element_t element;
          wrapped.calc_element (i, j, k, element);
          return element.get_corners ();};

      //! return coords of block vertexes by n_block index
      grd_ecl::fpoint3d_vector calc_element (const i_type_t index) const
        {return wrapped.calc_element (index).get_corners ();};

      ///////////////////////
      // ACCESS ARRAYS
      ///////////////////////


      //! get const int_to_ext
      const sp_i_array_t get_int_to_ext() const
        {return wrapped.get_int_to_ext ();};

      //! get const ext_to_int
      const sp_i_array_t get_ext_to_int() const
        {return wrapped.get_ext_to_int ();};

      //! get mesh elements volumes
      const sp_fp_array_t get_volumes () const
        {return wrapped.get_volumes ();};


      //! return depths of cell centers (length n_active_elements)
      const sp_fp_array_t get_depths () const
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
      int build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn, sp_i_array_t boundary_array)
        {return wrapped.build_jacobian_and_flux_connections (jacobian, flux_conn, boundary_array);};

      boost::python::list calc_element_tops ()
      {return wrapped.calc_element_tops();};

	
	  boost::python::list calc_element_center ()
      {return wrapped.calc_element_center();};

    
	//! init coord & zcorn from (nx, ny, nz, dx, dy, dz)
	//! return: first -- coord, second -- zcorn
	static std::pair< sp_fp_storage_array_t, sp_fp_storage_array_t >
	gen_coord_zcorn(i_type_t nx, i_type_t ny, i_type_t nz, sp_fp_storage_array_t dx, sp_fp_storage_array_t dy, sp_fp_storage_array_t dz) {
		return wrapped_t::gen_coord_zcorn(nx, ny, nz, dx, dy, dz);
	}

	static std::pair< sp_fp_storage_array_t, sp_fp_storage_array_t >
	refine_mesh(i_type_t& nx, i_type_t& ny, sp_fp_storage_array_t coord, sp_fp_storage_array_t zcorn, sp_fp_storage_array_t points) {
		return wrapped_t::refine_mesh(nx, ny, coord, zcorn, points);
	}

    ////////////////////
    // wrapped class
    ///////////////////

    private:

		typedef mesh_grdecl< strategy_t > wrapped_t;
		wrapped_t wrapped;
    };

};//namespace blue_sky
#endif // BS_MESH_GRDECL_H
