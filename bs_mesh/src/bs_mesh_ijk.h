#ifndef BS_MESH_IJK_H
#define BS_MESH_IJK_H
/*!
	\file bs_mesh_ijk.h
  \brief This file declares bs wrapper over mesh_ijk
  \author Mark Khait
	\date 2008-05-20
 */

#include "mesh_ijk.h"
#include "rs_smesh_iface.h"
#include "bs_flux_connections.h"

namespace blue_sky
  {

  
  class BS_API_PLUGIN bs_mesh_ijk : public rs_smesh_iface
    {

    ///////////////////////////////
    //  INTERNAL TYPE DECLARATION
    ///////////////////////////////

    public:
      ///////////////////////
      // BASE TYPES
      ///////////////////////
      typedef rs_smesh_iface                  base_t;

      ///////////////////////
      // OWN TYPES
      ///////////////////////
      typedef boost::array <grd_ecl::fpoint3d, 8>         fpoint3d_vector;

    public:
      //! blue-sky class declaration
      BLUE_SKY_TYPE_DECL(bs_mesh_ijk);

      //! default destructor
      ~bs_mesh_ijk ()	{};

      ///////////////////////
      // INITIALIZATION
      ///////////////////////


      //! init mesh properties
      void init_props (const sp_hdm_t hdm)
        {wrapped.init_props (hdm);};

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
      t_long get_n_active_elements () const
      {return wrapped.get_n_active_elements ();};

      //! return number of mesh elements
      t_long get_n_elements () const
      {return wrapped.get_n_elements ();};

      //! return number of mesh elements connections
      t_long get_n_connections () const
      {return wrapped.get_n_connections ();};

      //! return mesh dimensions range
      void get_dimensions_range (t_double &dim1_max, t_double &dim1_min,
                                         t_double &dim2_max, t_double &dim2_min,
                                         t_double &dim3_max, t_double &dim3_min) const
      {return wrapped.get_min_max_xyz (dim1_max, dim1_min, dim2_max, dim2_min, dim3_max, dim3_min);};

      //! get mesh dimensions
      base_t::index_point3d_t get_dimens ()
      {return wrapped.get_dimens();};

      //! return element size in all 3 dimensions
      void get_element_size (const t_long n_element, t_double &d_dim1, t_double &d_dim2, t_double &d_dim3) const
      {wrapped.get_block_dx_dy_dz(n_element, d_dim1, d_dim2, d_dim3);};

      //! return element size in 3rd dimension
      t_double get_element_dim3_size (const t_long n_element) const
      {return wrapped.get_block_dz(n_element);};

      //! return center point of an element by I, J and K mesh coordinates
      point3d_t get_element_center (const t_long i, const t_long j, const t_long k)const
      {return wrapped.get_center(i, j, k);};

      //! return center point of an element by internal element number
      point3d_t get_element_center (const t_long n_element)const
      {return wrapped.get_center(n_element);};

      //! return depth of mesh element
      t_double get_element_depth(const t_long n_element) const
      {return wrapped.get_depth(n_element);};

      //! return top of mesh element
      t_double get_element_dtop(const t_long n_element) const
      {return wrapped.get_dtop(n_element);};

      //! get element internal number by external
      t_long convert_ext_to_int (const t_long n_element) const
      {return wrapped.convert_ext_to_int(n_element);};

      //! get element external number by internal
      t_long get_element_int_to_ext (const t_long n_element) const
      {return wrapped.get_element_int_to_ext(n_element);};

      //! return I, J and K structured mesh coordinates of an element by internal number
      void get_element_int_to_ijk (const t_long n_element, t_long &i, t_long &j, t_long &k) const
      {wrapped.inside_to_XYZ(n_element, i, j, k);};

      //! return internal number of an element by I, J and K structured mesh coordinates
      t_long get_element_ijk_to_int (const t_long i, const t_long j, const t_long k) const
      {return wrapped.XYZ_to_inside(i, j, k);};

      //! return coords of block vertexes by IJK indexes
      grd_ecl::fpoint3d_vector calc_element (const t_long i, const t_long j, const t_long k) const
        {return wrapped.calc_element (i, j, k);};

      //! return coords of block vertexes by n_block index
      grd_ecl::fpoint3d_vector calc_element (const t_long index) const
        {return wrapped.calc_element (index);};

      ///////////////////////
      // ACCESS ARRAYS
      ///////////////////////


      //! get const int_to_ext
      const spv_long get_int_to_ext() const
      {return wrapped.get_int_to_ext ();};

      //! get const ext_to_int
      const spv_long get_ext_to_int() const
      {return wrapped.get_ext_to_int ();};

      //! get mesh elements volumes
      const spv_float get_volumes () const
      {return wrapped.get_volumes ();};


      //! return depths of cell centers (length n_active_elements)
      const spv_double get_depths () const
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
      int build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn, spv_long boundary_array)
      {return wrapped.build_jacobian_and_flux_connections (jacobian, flux_conn, boundary_array);};

	  spv_double
	  get_element_sizes (const t_long n_element) const
	  {
		double dx, dy, dz;
		get_element_size(n_element, dx, dy, dz);
		spv_double sizes = BS_KERNEL.create_object(v_double::bs_type());
		sizes->resize(3);
		(*sizes)[0] = dx;
		(*sizes)[1] = dy;
		(*sizes)[2] = dz;
		return sizes;
	  }
    ////////////////////
    // wrapped class
    ///////////////////

    private:

      mesh_ijk wrapped;
    };

};//namespace blue_sky
#endif // BS_MESH_IJK_H
