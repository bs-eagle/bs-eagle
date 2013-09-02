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

#define DEF_CELL_MERGE_THRESHOLD 0.8
#define DEF_BAND_THRESHOLD 0.2

namespace blue_sky
  {

  
  class BS_API_PLUGIN bs_mesh_grdecl : public rs_smesh_iface
    {

    ///////////////////////////////
    //  INTERNAL TYPE DECLARATION
    ///////////////////////////////
    public:
      ///////////////////////
      // BASE TYPES
      ///////////////////////
      typedef rs_smesh_iface                  base_t;
      typedef mesh_grdecl wrapped_t;

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
      void init_props (const sp_hdm_t hdm)
        {wrapped.init_props (hdm);};

      //! initialize int_to_ext indexation
      int init_int_to_ext()
        {return wrapped.init_int_to_ext ();};

      //! initialize int_to_ext indexation
      int init_ext_to_int()
        {return wrapped.init_ext_to_int();};

      //! init COORD & ZCORN via gen_coord_zcorn
      void init_props(t_long nx, t_long ny, t_long nz, spv_float dx, spv_float dy, spv_float dz) {
	return wrapped.init_props(nx, ny, nz, dx, dy, dz);
      }

      //! init COORD & ZCORN directly
      void init_props(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
	return wrapped.init_props(nx, ny, coord, zcorn);
      }


      //! reset internal state and release memory
      void clear() {
	wrapped.clear();
      }

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
      void get_dimensions_range (t_float& dim1_max, t_float &dim1_min,
        t_float &dim2_max, t_float &dim2_min,
        t_float &dim3_max, t_float &dim3_min) const
        {return wrapped.get_min_max_xyz (dim1_max, dim1_min, dim2_max, dim2_min, dim3_max, dim3_min);};

      //! get mesh dimensions
      index_point3d_t get_dimens ()
      {return wrapped.get_dimens();};

      //! return element size in all 3 dimensions
      void get_element_size (const t_long n_element, t_double &d_dim1, t_double &d_dim2, t_double &d_dim3) const
        {wrapped.get_block_dx_dy_dz(n_element, d_dim1, d_dim2, d_dim3);};
        
      spv_double get_element_sizes(const t_long n_element) const
        {return wrapped.get_element_sizes(n_element);};  

      //! return element size in 3rd dimension
      t_double get_element_dim3_size (const t_long n_element) const
        {return wrapped.get_block_dz(n_element);};
        
      //! return element size in 3rd dimension
      t_double get_element_dim3_size_ext (const t_long i, const t_long j, const t_long k) const
        {return wrapped.get_block_dz_ext(i, j, k);};

      //! return center point of an element
      point3d_t get_element_center (const t_long i, const t_long j, const t_long k)const
        {return wrapped.get_center(i, j, k);};

      //! return center point of an element by internal element number
      point3d_t get_element_center (const t_long n_element)const
        {return wrapped.get_center(n_element);};

      //! return depth of mesh element
      t_double get_element_depth(const t_long n_element) const
        {return wrapped.get_depth(n_element);};

      //! return number of mesh elements connections
      t_float get_element_dtop(const t_long n_element) const
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
        {  mesh_grdecl::element_t element;
          wrapped.calc_element (i, j, k, element);
          return element.get_corners ();};

      //! return coords of block vertexes by n_block index
      grd_ecl::fpoint3d_vector calc_element (const t_long index) const
        {return wrapped.calc_element (index).get_corners ();};

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

      //! get cell volumes
      spv_float get_cell_volumes(const t_long Nx, const t_long Ny, const t_long Nz) const {
          return wrapped.get_cell_volumes(Nx, Ny, Nz);
      }

      //! return depths of cell centers (length n_active_elements)
      const spv_float get_depths () const
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

      //! find well`s trajectories and mesh cells intersection
      int intersect_trajectories (sp_well_pool_t well_pool) {return wrapped.intersect_trajectories(well_pool);};

      boost::python::list calc_element_tops ()
      {return wrapped.calc_element_tops();};

      // same as calc_element_tops, but only return tops coordinates
      spv_float calc_cells_vertices() {
		  return wrapped.calc_cells_vertices();
	  }
      spv_float calc_cells_vertices_xyz() {
		  return wrapped.calc_cells_vertices_xyz();
	  }

      boost::python::list calc_element_center ()
      {return wrapped.calc_element_center();};

    
	//! init coord & zcorn from (nx, ny, nz, dx, dy, dz)
	//! return: first -- coord, second -- zcorn
	static std::pair< spv_float, spv_float >
	gen_coord_zcorn(t_long nx, t_long ny, t_long nz, spv_float dx, spv_float dy, spv_float dz,
			t_float x0 = 0, t_float y0 = 0, t_float z0 = 0)
	{
		return wrapped_t::gen_coord_zcorn(nx, ny, nz, dx, dy, dz, x0, y0, z0);
	}

  //! init coord
	static spv_float
	gen_coord2(spv_float x, spv_float y)
	{
		return wrapped_t::gen_coord2(x, y);
	}

	static std::pair< spv_float, spv_float >
	refine_mesh_deltas(t_long& nx, t_long& ny, spv_float coord, spv_float points,
			spv_long hit_idx = NULL,
			t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD)
	{
		return wrapped_t::refine_mesh_deltas(nx, ny, coord, points, hit_idx, cell_merge_thresh, band_thresh);
	}

	static std::pair< spv_float, spv_float >
	refine_mesh_deltas(t_long& nx, t_long& ny, spv_float coord,
			spv_long points_pos, spv_float points_param,
			spv_long hit_idx = NULL,
			t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD)
	{
		return wrapped_t::refine_mesh_deltas(nx, ny, coord, points_pos, points_param,
				hit_idx, cell_merge_thresh, band_thresh);
	}

	static std::pair< spv_float, spv_float >
	refine_mesh(t_long& nx, t_long& ny, spv_float coord, spv_float zcorn, spv_float points,
			spv_long hit_idx = NULL,
			t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD)
	{
		return wrapped_t::refine_mesh(nx, ny, coord, zcorn, points, hit_idx, cell_merge_thresh, band_thresh);
	}

	static std::pair< spv_float, spv_float >
	refine_mesh(t_long& nx, t_long& ny, spv_float coord, spv_float zcorn,
			spv_long points_pos, spv_float points_param,
			spv_long hit_idx = NULL,
			t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD)
	{
		return wrapped_t::refine_mesh(nx, ny, coord, zcorn, points_pos, points_param,
				hit_idx, cell_merge_thresh, band_thresh);
	}

	static std::pair< sp_fp_storage_array_t, sp_fp_storage_array_t >
	refine_mesh(i_type_t& nx, i_type_t& ny, sp_fp_storage_array_t coord, sp_fp_storage_array_t zcorn, sp_fp_storage_array_t points,
			sp_i_array_t hit_idx = NULL,
			fp_type_t cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, fp_type_t band_thresh = DEF_BAND_THRESHOLD)
	{
		return wrapped_t::refine_mesh(nx, ny, coord, zcorn, points, hit_idx, cell_merge_thresh, band_thresh);
	}

    ////////////////////
    // wrapped class
    ///////////////////

    // direct access to wrapped class
    const wrapped_t& get_wrapped() const {
      return wrapped;
    }

    wrapped_t& get_wrapped() {
      return wrapped;
    }

    private:

		wrapped_t wrapped;
    };

};//namespace blue_sky
#endif // BS_MESH_GRDECL_H
