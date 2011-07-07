#ifndef MESH_GRDECL_H
#define MESH_GRDECL_H
/**
	\file mesh_grdecl.h
 	\brief This file declare class for working with grid_ecllipse
	\author Iskhakov Ruslan
	\date 2008-05-20 */

#include "rs_smesh_base.h"
#include "flux_connections_iface.h"

#include "mesh_element3d.h"

#ifdef _HDF5_MY //!< using HDF5 or not
#include "H5Cpp.h"
#endif	//_HDF5_MY

#ifndef EPS_SQ
#define EPS_SQ 1.0e-5 //!< tolerance for calculating of area
#endif // EPS_SQ

#ifndef HDF5_MAX_SPACE
#define HDF5_MAX_SPACE 40 //!< maxinum size of space that we can read from ASCII file
#endif

#define DEF_CELL_MERGE_THRESHOLD 0.8
#define DEF_BAND_THRESHOLD 0.2

//! class for basic work with mesh based on ZCORN&COORD and tpfa calculating

class BS_API_PLUGIN mesh_grdecl : public  rs_smesh_base
  {


    friend struct build_jacobian_and_flux;
    template <class L> friend struct build_jacobian_rows_class;
    template <class L> friend struct build_jacobian_cols_class;
  private:
	  struct inner;
	  st_smart_ptr< inner > pinner_;
  
//+++++++++++++++++++++++++++++++++++++++++++
//  INTERNAL TYPE DECLARATION
//===========================================
  public:
    ///////////////////////
    // BASE TYPES
    ///////////////////////
    typedef rs_smesh_base                   base_t;


    ///////////////////////
    // OWN TYPES
    ///////////////////////

    typedef mesh_element3d                  element_t;
    typedef smart_ptr <element_t, true>                 sp_element_t;
    typedef element_t::corners_t               corners_t;
    typedef element_t::plane_t                 plane_t;

    typedef grd_ecl::fpoint3d                  fpoint3d_t;
    typedef grd_ecl::fpoint2d                  fpoint2d_t;
    typedef grd_ecl::quadrangle_t              quadrangle_t;
    
    typedef boost::array <t_long, 8>         element_zcorn_t_long;
    typedef boost::array <t_long, 4>         plane_zcorn_t_long;
    typedef boost::array <t_double, 3>                point3d_t;

    typedef blue_sky::smart_ptr <blue_sky::FRead, true> sp_fread_t;

//-------------------------------------------
//  METHODS
//===========================================
  public:
    //! default constructor
    mesh_grdecl ();

    //! default destructor
    ~mesh_grdecl () {};

    //! read zcorn from FRead
    //void read_zcorn(const sp_fread_t &r);
    
    //! init arrays of properties
    void init_props (const sp_hdm_t hdm);

	//! init COORD & ZCORN via gen_coord_zcorn
	void init_props(t_long nx, t_long ny, t_long nz, spv_float dx, spv_float dy, spv_float dz);

	//! init COORD & ZCORN directly
	void init_props(t_long nx, t_long ny, spv_float coord, spv_float zcorn);

	//! init COORD & ZCORN from (nx, ny, nz, dx, dy, dz)
	//! (x0, y0, z0) - beginning of coordinate system (mesh offset)
	//! return: first -- coord, second -- zcorn
	static std::pair< spv_float, spv_float >
	gen_coord_zcorn(t_long nx, t_long ny, t_long nz, spv_float dx, spv_float dy, spv_float dz,
		t_float x0 = 0, t_float y0 = 0, t_float z0 = 0);

	//! return refined dx and dy for given COORD
	static std::pair< spv_float, spv_float >
	refine_mesh_deltas(t_long& nx, t_long& ny, spv_float coord, spv_float points,
		spv_long hit_idx = NULL,
		t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD);

	//! return refined dx and dy for given COORD
	//! points are given in (i,j) format
	static std::pair< spv_float, spv_float >
	refine_mesh_deltas(t_long& nx, t_long& ny, spv_float coord,
		spv_long points_pos, spv_float points_param,
		spv_long hit_idx = NULL,
		t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD);

	//! refine_mesh = refine_mesh_deltas + gen_coord_zcorn
	//! return COORD & ZCORN of refined mesh
	static std::pair< spv_float, spv_float >
	refine_mesh(t_long& nx, t_long& ny, spv_float coord, spv_float zcorn, spv_float points,
		spv_long hit_idx = NULL,
		t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD);

	//! refine_mesh = refine_mesh_deltas + gen_coord_zcorn
	//! points are given in (i,j) format
	//! return COORD & ZCORN of refined mesh
	static std::pair< spv_float, spv_float >
	refine_mesh(t_long& nx, t_long& ny, spv_float coord, spv_float zcorn,
		spv_long points_pos, spv_float points_param,
		spv_long hit_idx = NULL,
		t_double cell_merge_thresh = DEF_CELL_MERGE_THRESHOLD, t_double band_thresh = DEF_BAND_THRESHOLD);

    //! get vertex of cube [i,j,k]
    void calc_element (const t_long i, const t_long j, const t_long k, element_t &element) const;
    
    //! get vertex of cube [index]
    element_t calc_element (const t_long index) const;

    //! get coords && zcorn from file_name
    bool file_open_cube(const char* file_name);

    //! get actnum array
    bool file_open_actnum(const char* file_name);

    /**
    * \brief Reads "coord" and "zcorn" arrays from ascii file
    * \param filename -- path to the ascii file
    * \return true if success, false otherwise
    */
#ifdef _HDF5_MY
    bool file_open_cube_with_hdf5_swap(const char* file_name);
    /**
    * \brief Create 1-dimensional t_double dataset in hdf5 file
    * \param arr - dataset name
    * \param file_hdf5 - hdf5 file
    * \param dataset
    * \return 0 if success
    */
    int create_array_hdf5(const char *dataset_name, H5::H5File &file_hdf5, H5::DataSet **dataset);
    /**
    * \brief Adds specified array to the specified hdf5 dataset
    * \param arr - array to add
    * \param arr_length - array length
    * \param dataset - hdf5 dataset
    * \return 0 if success
    */
    int append_array_hdf5(const t_float *arr, size_t arr_length, H5::DataSet *dataset);
    bool file_open_activs_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt);
    bool file_open_cube_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt);
#endif
    ///////////////


    //-------------------------------------------
    //  INHERITED FUNCTIONS
    //===========================================
    //! make indexation - create proxy and non-proxy array using info of minpv
    //! return number of non active blocks (less than mpv)
    int init_ext_to_int();


    /*!	\brief allocate jacobian and fill structure
    	\param n_block_size matrix block size
    	\param jacobian jacobian (like adjacency matrix)
    	\param flux_conn connection with calculated transmissibility
    	\return 0 if success

    	we optimize this function by skipping butting cells and change bypass method (Z->Y->X => X->Y->Z);
    	in all cases we using local indexing.*/
    int build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn, 
                                             spv_long boundary_array);

    /*!	\brief allocate jacobian and fill structure
    	\param n_block_size matrix block size
    	\param jacobian jacobian (like adjacency matrix)
    	\param flux_conn connection with calculated transmissibility
    	\param boundary_array array of block, which obey boundary condition
    	\return number of boundary cell

    	we optimize this function by skipping butting cells and change bypass method (Z->Y->X => X->Y->Z);
    	boundary condition = active block, which have non-active neighbors;
    	in all cases we using local indexing.*/
    int build_jacobian_and_flux_connections_add_boundary (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn,
                                                          spv_long boundary_array);

	boost::python::list calc_element_tops ();

	boost::python::list calc_element_center ();

	// same as calc_element_tops, but only return tops coordinates
	spv_float calc_cells_vertices();

    /*!	\brief  find neighbours (adjacency matrix)
    	\param neig_matrix - bcsr adjacency matrix
    	\return 0 if success */
    int find_neighbours(sp_bcsr_t /*neig_matrix*/) {return 0;};

    /*! \brief get_block_dx_dy_dz
        \param n_elem - block number
        \param (return value) dx, dy, dz - average size of current block

      average_size dx = fabs[ average(right_side by x) - average (left_side by x)]
      average_size dy = fabs[ average(top_side by y) - average (bottom_side by y)]
      average_size dz = fabs[ average(upper_side by z) - average (lower_side by z)]*/
    void get_block_dx_dy_dz(t_long n_elem, t_double &dx, t_double &dy, t_double &dz) const;

    //! short  analog for previous function
    t_double get_block_dx(t_long n_elem) const;

    t_double get_block_dy(t_long n_elem) const;

    t_double get_block_dz(t_long n_elem) const;
    
    t_double get_depth(t_long n_elem) const
      {
        return (*depths)[n_elem];
      }; 

    t_double get_dtop(t_long n_elem) const;

    //! function for filling net by testing data
    void generate_array();
    
    //! check mesh data
    void check_data () const;

  protected:
    /*! \brief fill given array with block centers (only activ block) and according local bypass
    	\return 0 if success */
    int calc_depths ();

    //! define type of side for block [i,j,k] and his neighbors [i1,j1,k1]
    int define_side(const t_long i, const t_long j,const  t_long k, const t_long i1, const t_long j1, const t_long k1) const;

    //! crossing coord&z - real coordinat of one point
    void calc_corner_point (const t_float z, const t_float *coord, fpoint3d_t &p) const;

    /*!	\brief is block[i,j,k] have been already splicing
    	\return true if already splicing, else - false */
    bool is_small(const t_long i, const t_long j, const t_long k, t_double eps) const;

    /*!	\brief calc pore volume for block[i,j,k]
        \return pore volume */
    t_double calc_block_volume(const t_long i, const t_long j, const t_long k) const ;

    /*! \brief splice 2 blocks (block [i,j,k] absorbs block [i,j,k1]) recalculating zcorn, perms, pore volume, NTG*/
    void splice_two_blocks (const t_long i, const t_long j, const t_long k, const t_long k1) ;

    /*! \brief check if 2 blocks (block [i,j,k] and [i,j,k1]) are close enough to be coupled*/
    bool are_two_blocks_close (const t_long i, const t_long j, const t_long k, const t_long k1) ;
    
    /*! \brief check if every active cell in mesh is adjacent*/
    bool check_adjacency(int shift_zcorn = 0);

  public:
    //! calculate projection of area for all axis (for transmissibility and others calculation)
    void get_plane_crossing_projection_on_all_axis(const plane_t &side1, const plane_t &side2, fpoint3d_t &A)const;

    plane_zcorn_t_long get_plane_zcorn_index (element_zcorn_t_long &element, element_plane_orientation_t orientation) const;

    void get_element_zcorn_index (const t_long i, const t_long j, const t_long k, element_zcorn_t_long& element) const;


  protected:


    //! get arithmetic center of cube, which elements are indexes from zcorn_array[]
    t_double get_center_zcorn(const element_zcorn_t_long &element) const;

    /*! \brief splicing block if it's too small
        \param volumes_temp - array of volumes
    	  \return number of splicing blocks

    	if block too small, we change it's params (zcorn) => it has volume 0, but neighbor block will be grow;
    	also this too small block will be non-active (change actnum_array)*/
    int  splicing(stdv_float & volumes_temp);

  public:

    /*! \brief calculate transmissbility between two blocks ()
    	\param ext_index1 - global indexation (X->Y->Z) for first block
    	\param ext_index2 - global indexation (X->Y->Z) for second block
    	\param plane1 - contact plane of block1
    	\param center1 - center of block1
    	\param center2 - center of block2
    	\param d_dir - direction of transmissibility
    	\param plane2 - contact plane of block2 (if 0, then blocks are fully adjacent and contact area calculation is simple)
    	\return transmissiblity value between 2 blocks */
      t_double calc_tran (const t_long ext_index1, const t_long ext_index2, const plane_t &plane1, 
                          const fpoint3d_t &center1, const fpoint3d_t &center2, direction d_dir, plane_t *plane2  = 0) const;                          


//-------------------------------------------
//  VARIABLES
//===========================================
  public:

    // TODO: replace with seq_vector
  public:
    t_float *coord_array;	          //!< COORD array
    t_float *zcorn_array;				    //!< ZCORN array
  };


#endif //MESH_GRDECL_H
