#ifndef MESH_GRDECL_H
#define MESH_GRDECL_H
/**
	\file mesh_grdecl.h
 	\brief This file declare class for working with grid_ecllipse
	\author Iskhakov Ruslan
	\date 2008-05-20 */

#include "rs_smesh_base.h"
#include "flux_connections_iface.h"

#include "fpoint3d.h"

#ifdef _HDF5_MY //!< using HDF5 or not
#include "H5Cpp.h"
#endif	//_HDF5_MY

#ifndef EPS_SQ
#define EPS_SQ 1.0e-5 //!< tolerance for calculating of area
#endif // EPS_SQ

#ifndef HDF5_MAX_SPACE
#define HDF5_MAX_SPACE 40 //!< maxinum size of space that we can read from ASCII file
#endif


namespace rs_mesh_detail {

  //! enum for define type of crossing between planes of cells
  enum cross_section
  {
    empty, //!< no intersection
    left, right,
    top, bottom,
    upper, lower,
  };

  //! element of coord line (2 point - start and begin)
  struct coordElem
  {
    grd_ecl::fpoint3d pStart;	//!< start point of coord
    grd_ecl::fpoint3d pEnd;		//!< end point of coord

    //!< default constructor
    coordElem() 
    {
    }	

    //! constructor
    coordElem(const grd_ecl::fpoint3d &apStart, const grd_ecl::fpoint3d &apEnd)
    {
      pStart = apStart;
      pEnd = apEnd;
    }
  };

} // namespace mesh_detail


//! class for basic work with mesh based on ZCORN&COORD and tpfa calculating
template<class strategy_t>
class BS_API_PLUGIN mesh_grdecl : public  mesh_rs<strategy_t>
  {
  
//+++++++++++++++++++++++++++++++++++++++++++
//  INTERNAL TYPE DECLARATION
//===========================================
  public:
    ///////////////////////
    // BASE TYPES
    ///////////////////////
    typedef mesh_rs <strategy_t>                        base_t;

    typedef typename base_t::index_t                    index_t;
    typedef typename base_t::item_t                     item_t;

    typedef typename base_t::index_array_t              index_array_t;
    typedef typename base_t::item_array_t               item_array_t;

    typedef typename base_t::sp_bcsr_t                  sp_bcsr_t;
    typedef typename base_t::sp_idata_t                 sp_idata_t;
    typedef typename base_t::sp_flux_conn_iface_t       sp_flux_conn_iface_t;
    
    //typedef typename base_t::i_map_t                    i_map_t;
    //typedef typename base_t::d_map_t                    d_map_t;

    ///////////////////////
    // OWN TYPES
    ///////////////////////

    typedef grd_ecl::fpoint3d                           fpoint3d;
    typedef std::vector<grd_ecl::fpoint2d>              g_fpoint2d_vector;
    typedef typename strategy_t::rhs_item_array_t       rhs_item_array_t;
    typedef typename array_float16_t::value_type        pool_item_t;

    typedef grd_ecl::fpoint3d_vector                    ijk_cube_t;
    typedef boost::array <index_t, 4>                   side_t;
    typedef boost::array <fpoint3d, 4>                  point_side_t;
    typedef boost::array <item_t, 3>                    center_t;
    typedef boost::array <index_t, 8>                   cube_index_t;

    typedef blue_sky::smart_ptr <blue_sky::FRead, true> sp_fread_t;


//-------------------------------------------
//  METHODS
//===========================================
  public:
    //! default constructor
    mesh_grdecl () {};

    //! default destructor
    ~mesh_grdecl () {};

    //! read zcorn from FRead
    //void read_zcorn(const sp_fread_t &r);
    
    //! init arrays of properties
    void init_props (const sp_idata_t &idata);
    
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
    * \brief Create 1-dimensional item_t dataset in hdf5 file
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
    int append_array_hdf5(const item_t *arr, size_t arr_length, H5::DataSet *dataset);
    bool file_open_activs_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt);
    bool file_open_cube_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt);
#endif
    ///////////////

    center_t
    get_center (index_t i, index_t j, index_t k) const;

    center_t
    get_center (index_t n_block) const;

    //! get volume of current fpoint3d cube
    float get_volume_cube (const grd_ecl::fpoint3d_vector &cube) const;

    //! get vertex of cube (i,j,k) and center
    //! length of (fpoint3d_vector) = 8
    grd_ecl::fpoint3d_vector top_cube (const index_t i, const index_t j, const index_t k) const;

    //! get vertex of cube (i,j,k) and center
    //! length of (fpoint3d_vector) = 8
    grd_ecl::fpoint3d_vector top_cube (const index_t index) const;

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
    int build_jacobian_and_flux_connections (const sp_bcsr_t &jacobian, const sp_flux_conn_iface_t &flux_conn, 
                                             index_array_t &boundary_array);

    /*!	\brief allocate jacobian and fill structure
    	\param n_block_size matrix block size
    	\param jacobian jacobian (like adjacency matrix)
    	\param flux_conn connection with calculated transmissibility
    	\param boundary_array array of block, which obey boundary condition
    	\return number of boundary cell

    	we optimize this function by skipping butting cells and change bypass method (Z->Y->X => X->Y->Z);
    	boundary condition = active block, which have non-active neighbors;
    	in all cases we using local indexing.*/
    int build_jacobian_and_flux_connections_add_boundary (const sp_bcsr_t &jacobian, const sp_flux_conn_iface_t &flux_conn,
                                                          index_array_t &boundary_array);

    /*!	\brief  find neighbours (adjacency matrix)
    	\param neig_matrix - bcsr adjacency matrix
    	\return 0 if success */
    int find_neighbours(sp_bcsr_t &neig_matrix);

    /*! \brief get_block_dx_dy_dz
        \param n_elem - block number
        \param (return value) dx, dy, dz - average size of current block

      average_size dx = fabs[ average(right_side by x) - average (left_side by x)]
      average_size dy = fabs[ average(top_side by y) - average (bottom_side by y)]
      average_size dz = fabs[ average(upper_side by z) - average (lower_side by z)]*/
    void get_block_dx_dy_dz(index_t n_elem, item_t &dx, item_t &dy, item_t &dz) const;

    //! short  analog for previous function
    float get_block_dx(index_t n_elem) const;

    float get_block_dy(index_t n_elem) const;

    float get_block_dz(index_t n_elem) const;
    
    float get_depth(index_t n_elem) const
      {
        return depths[n_elem];
      }; 

    float get_dtop(index_t n_elem) const;

    //! function for filling net by testing data
    void generate_array();
    
    //! check mesh data
    void check_data () const;

  protected:
    /*! \brief fill given array with block centers (only activ block) and according local bypass
    	\return 0 if success */
    int calc_depths ();

    //! define type of side for block [i,j,k] and his neighbors [i1,j1,k1]
    int define_side(const index_t i, const index_t j,const  index_t k, const index_t i1, const index_t j1, const index_t k1) const;

    //! crossing coord&z - real coordinat of one point
    fpoint3d cross_coord(const item_t z, const pool_item_t *coord) const;

    /*!	\brief is block[i,j,k] have been already splicing
    	\return true if already splicing, else - false */
    bool is_small(const index_t i, const index_t j, const index_t k, item_t eps) const;

    /*!	\brief calc pore volume for block[i,j,k]
        \return pore volume */
    item_t calc_block_volume(const index_t i, const index_t j, const index_t k) const ;

    /*! \brief splice 2 blocks (block [i,j,k] absorbs block [i,j,k1]) recalculating zcorn, perms, pore volume, NTG*/
    void splice_two_blocks (const index_t i, const index_t j, const index_t k, const index_t k1) ;

    /*! \brief check if 2 blocks (block [i,j,k] and [i,j,k1]) are close enough to be coupled*/
    bool are_two_blocks_close (const index_t i, const index_t j, const index_t k, const index_t k1) ;

  public:
    //! calculate projection of area for all axis (for transmissibility and others calculation)
    fpoint3d get_side_crossing_projection_on_all_axis(const point_side_t &side1, const point_side_t &side2)const;

  protected:

    ///*! \brief get elems from cube which from selected side
    //	\param i_side side, which we take from cube
    //	\param cube source vector of vertex (or indexex)
    //	\param is_anti if we take i_side, is_anti must be false; if we take antiside (left->right) is_anti must be true
    //	\return selected side, defined by params */
    //template <class array_t>
    //array_t get_side(cross_section i_side, const array_t &cube, const bool is_anti)const;

    //! get arithmetic center of cube, which elements are indexes from zcorn_array[]
    item_t get_center_zcorn(const cube_index_t &cube) const;

    /*! \brief splicing block if it's too small
        \param volumes_temp - array of volumes
    	  \return number of splicing blocks

    	if block too small, we change it's params (zcorn) => it has volume 0, but neighbor block will be grow;
    	also this too small block will be non-active (change actnum_array)*/
    int  splicing(item_array_t& volumes_temp);

    /*! \brief add one element[i,j,k] in rows_ptr*/
    void change_row(const index_t i, const index_t j, const index_t k, index_array_t &rows_ptr);

    /*! \brief add index2 into cols_ind by address curIndex[index1]*/
    void change_col(const index_t index1, const index_t index2, index_array_t &cols_ind, index_array_t& curIndex);

    /*! \brief add index2 into cols_ind by address curIndex[index1] and change m&p_memory block indexation
      \param is_m_memory define what array use to push_back index2
      \param connection_number define what number of connection we put in*/
    void change_col_add(const index_t index1, const index_t index2, index_array_t &cols_ind,
                        index_array_t& curIndex, index_array_t& m_memory, index_array_t& p_memory, bool is_m_memory,
                        bool is_need_to_add, const index_t connection_number, const index_array_t &rows_ptr);

    /*! \brief set pointers to index1 and index2 connection in jacobian
      \param connection_number define what number of connection we put in*/
    void set_plus_minus (const index_t index1, const index_t index2, index_array_t &cols_ind, index_array_t &m_memory, index_array_t &p_memory,
                         const index_t connection_number_minus, const index_t connection_number_plus,  const index_array_t &rows_ptr);

  public:
    /*!	\brief is 2 side crossing or not (using zcorn dependences between blocks)
    	\param side1 - fpoint3d vector of selected side (block 1)
    	\param side2 - fpoint3d vector of selected side (block 2)
    	\param x_dir - true if dir = {left, right}, else (dir = {top, bottom}) = false
    	\return true if 2 side have common intersection area*/
    bool is_2side_crossing(const side_t &side1, const side_t &side2, bool x_dir);

    /*! \brief calculate transmissbility between two blocks
    	\param index1 - global indexation (X->Y->Z) for first block
    	\param inded2 - global indexation (X->Y->Z) for second block
    	\param A - intersection area projection on all axis
    	\param D1 - difference between center of side1 && center of block1
    	\param side2 - fpoint3d vector of vertex (selected side)
    	\param center2 - center of block2
    	\param d_dir - direction of transmissibility
    	\return transmissiblity value between 2 blocks */
    item_t calculate_tran(const index_t index1, const index_t index2, const fpoint3d& A, const fpoint3d& D1, const point_side_t &side2,const fpoint3d &center2, direction d_dir) const;


//-------------------------------------------
//  VARIABLES
//===========================================
  public:
    using base_t::init_int_to_ext;
    using base_t::get_n_active_elements;

    using base_t::nx;
    using base_t::ny;
    using base_t::nz;
    using base_t::minpv;
    using base_t::minsv;
    using base_t::max_thickness;
    using base_t::min_x;
    using base_t::min_y;
    using base_t::min_z;
    using base_t::max_x;
    using base_t::max_y;
    using base_t::max_z;
    using base_t::depths;
    using base_t::ext_to_int;
    using base_t::int_to_ext;
    using base_t::sp_actnum;
    using base_t::sp_permx;
    using base_t::sp_permy;
    using base_t::sp_permz;
    using base_t::sp_poro;
    using base_t::sp_ntg;
    using base_t::sp_multx;
    using base_t::sp_multy;
    using base_t::sp_multz;
    using base_t::n_elements;    
    using base_t::n_active_elements;
    using base_t::n_connections;
    using base_t::darcy_constant;
    using base_t::volumes;


    // TODO: replace with seq_vector
  public:
    array_float16_t coord_array;	          //!< COORD array
    array_float16_t zcorn_array;				    //!< ZCORN array
  };

//! get arithmetic center of cube
template<class array_t>
typename array_t::value_type 
get_cube_center(const array_t &cube)
{
  typename array_t::value_type res;

  for (size_t i = 0; i < cube.size(); i++)
    {
      res += cube[i];
    }

  return res / cube.size();
}

//! get vertex's index in zcorn_array of block[i,j,k]
template<class index_t>
boost::array <index_t, 8>
top_cube_index(index_t i, index_t j, index_t k, index_t nx, index_t ny) 
{
  boost::array <index_t, 8> cubeIndex;
  //define index
  index_t index1 = i * 2 + j * 4 * nx + k * 8 * nx * ny;
  index_t index2 = index1 + 4 * nx * ny;

  cubeIndex[0] = index1;
  cubeIndex[1] = index1 + 1;
  cubeIndex[2] = index1 + 2 * nx;
  cubeIndex[3] = index1 + 2 * nx + 1;

  cubeIndex[4] = index2;
  cubeIndex[5] = index2 + 1;
  cubeIndex[6] = index2 + 2 * nx;
  cubeIndex[7] = index2 + 2 * nx + 1;

  return cubeIndex;
}

template <class array_t>
boost::array <typename array_t::value_type, 4> 
get_side(rs_mesh_detail::cross_section i_side, const array_t &cube, const bool is_anti)
{
  using namespace rs_mesh_detail;

  if (is_anti)
    {
      switch (i_side)
        {
        case left :
          i_side = right;
          break;
        case upper:
          i_side = lower;
          break;
        case top:
          i_side = bottom;
          break;
        case right:
          i_side = left;
          break;
        case bottom:
          i_side = top;
          break;
        case lower:
          i_side = upper;
          break;
        default:
          bs_throw_exception ("Invalid i_size value");
        };
    }

  boost::array <typename array_t::value_type, 4> side;
  switch (i_side)
    {
    case left :
      side[0] = cube[2];
      side[1] = cube[0];
      side[2] = cube[6];
      side[3] = cube[4];
      break;
    case right:
      side[0] = cube[3];
      side[1] = cube[1];
      side[2] = cube[7];
      side[3] = cube[5];
      break;
    case upper:
      side[0] = cube[6];
      side[1] = cube[4];
      side[2] = cube[7];
      side[3] = cube[5];
      break;
    case lower:
      side[0] = cube[2];
      side[1] = cube[0];
      side[2] = cube[3];
      side[3] = cube[1];
      break;
    case top:
      side[0] = cube[4];
      side[1] = cube[5];
      side[2] = cube[0];
      side[3] = cube[1];
      break;
    case bottom:
      side[0] = cube[6];
      side[1] = cube[7];
      side[2] = cube[2];
      side[3] = cube[3];
      break;
    default:
      bs_throw_exception ("Invalid i_size value");;
    }

  return side;
}

#endif //MESH_GRDECL_H
