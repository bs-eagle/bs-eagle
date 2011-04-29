#ifndef RS_SMESH_H
#define RS_SMESH_H
/*!
	\file rs_smesh_base.h
	\brief This file declare base class for reservoir simulation structured 3D meshes
	\author Iskhakov Ruslan
	\date 2008-05-20
 */

#include "flux_connections.h"
#include "rs_mesh_base.h"

using namespace blue_sky;

  //! enum for define direction
  enum direction
  {
    along_dim1 = 0,
    along_dim2,
    along_dim3
  };



class BS_API_PLUGIN rs_smesh_base : public rs_mesh_base
  {


  //-----------------------------------------
  // TYPES
  //-----------------------------------------

  public:
    ///////////////////////
    // BASE TYPES
    ///////////////////////
    typedef rs_mesh_base                    base_t;


    ///////////////////////
    // OWN TYPES
    ///////////////////////

    typedef rs_smesh_base                    this_t;
    typedef boost::array <t_long, 3>                   index_point3d_t;
    typedef std::pair<t_long, t_long>                elem_index;
    //typedef boost::array <grd_ecl::fpoint3d, 8>         fpoint3d_vector;

//-------------------------------------------
//  METHODS
//===========================================
  public:
    //! default constructor
    rs_smesh_base ();

    //! default destructor
    virtual ~rs_smesh_base ()	{};

    //! init mesh
    void init_props (const sp_hdm_t hdm);

    //! check mesh data
    void check_data() const;

    //! get min&max xyz coordinates
    void get_min_max_xyz (t_float &amax_x, t_float &amin_x, t_float &amax_y, t_float &amin_y, t_float &amax_z, t_float &amin_z) const
      {
        amax_x = max_x;
        amin_x = min_x;
        amax_y = max_y;
        amin_y = min_y;
        amax_z = max_z;
        amin_z = min_z;
      }

    index_point3d_t get_dimens () const
      {
        index_point3d_t dims;
        dims[0] = nx;
        dims[1] = ny;
        dims[2] = nz;
        return dims;
      }
    
    t_long get_nx () { return nx;};

    t_long get_ny () { return ny;};

    t_long get_nz () { return nz;};

    /*	\brief count how many activ blocks in layers, forming according direction
    	\param dir - direction of layers dividing
    	\param elem_in_layers - (return value) vector of active block's number in layer[]
    	\return number of layers
    	implementation non-dependece from type of mesh (just from actnum_array)		 */
    int get_elems_n_in_layers(const direction d_dir, stdv_int &elem_in_layers) const;

    ////////indexation functions

    //! return index of block[i,j,k] - same as XYZ_to_inside
    t_long get_n_block (const t_long i, const t_long j, const t_long k) const
      {
        return XYZ_to_inside (i, j, k);
      }

    //! convert global [i,j,k] to local indexing
    t_long inline XYZ_to_inside(const t_long i, const t_long j, const t_long k) const
      {
        if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz)
          return -1;
        return (*base_t::base_t::ext_to_int)[i + j*nx + k*nx*ny];
      }

    //! convert global [i+j*nx+k*nx*ny] to local indexing
    t_long XYZ_to_inside(const t_long i) const
      {
        return (*base_t::base_t::ext_to_int)[i];
      }

    //! convert local index to global [i,j,k]
    void inside_to_XYZ(const t_long index, t_long &i1, t_long &j1, t_long &k1)const;

    //! convert local index to global [i+j*nx+k*nx*ny]
    t_long inside_to_XYZ(const t_long index) const
      {
        return (*base_t::base_t::int_to_ext)[index];
      }

    /*
    //! make ext_to_int of num (int_to_ext)  (like original_elements_num in old mesh)
    int init_int_to_ext();
    */
    
    //-------------------------------------------
    //  VIRTUAL FUNCTIONS
    //===========================================
    
    /*
    virtual
    grd_ecl::fpoint3d_vector calc_element (const i_type_t i, const i_type_t j, const i_type_t k) const
      {
        BS_ASSERT (false && "PURE CALL");
        static grd_ecl::fpoint3d_vector dummy;
        return dummy;
      };
    */
    
    /*! \brief get_block_dx_dy_dz
        \param n_elem - block number
        \param (return value) dx, dy, dz - average size of current block

      average_size dx = fabs[ average(right_side by x) - average (left_side by x)]
      average_size dy = fabs[ average(top_side by y) - average (bottom_side by y)]
      average_size dz = fabs[ average(upper_side by z) - average (lower_side by z)]*/
    virtual void get_block_dx_dy_dz(t_long /*n_elem*/, t_double & /*dx*/, t_double & /*dy*/, t_double & /*dz*/) const
      {
        BS_ASSERT (false && "PURE CALL");
      };

    //! short  analog for previous function
    virtual t_double get_block_dx(t_long /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

    virtual t_double get_block_dy(t_long /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

    virtual t_double get_block_dz(t_long /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

    virtual t_double get_dtop(t_long /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      };

    // TODO: WTF??!!!
    virtual t_long get_coord_size() const
    {
      BS_ASSERT (false && "PURE CALL");
      return 0;
    }
    virtual t_long get_zcorn_size() const
    {
      BS_ASSERT (false && "PURE CALL");
      return 0;
    }

    virtual boost::array <t_double, 3>
    get_center (t_long /*i*/, t_long /*j*/, t_long /*k*/) const
      {
        BS_ASSERT (false && "PURE CALL");

        static boost::array <t_float, 3> dummy;
        return dummy;
      }

    virtual boost::array <t_double, 3>
    get_center (t_long /*n_element*/) const
      {
        BS_ASSERT (false && "PURE CALL");

        static boost::array <t_float, 3> dummy;
        return dummy;
      }

//-------------------------------------------
//  VARIABLES
//===========================================
  protected:
    t_long nx, ny, nz;	//!< grid dimens

    t_double min_x, max_x; //!< min&max x coordinates
    t_double min_y, max_y; //!< min&max y coordinates
    t_double min_z, max_z; //!< min&max z coordinates

    spv_float permx_array;	  //!< smart_ptr on permx_array array
    spv_float permy_array;	  //!< smart_ptr on permy_array array
    spv_float permz_array;	  //!< smart_ptr on permz_array array

    spv_float multx_array;	  //!< smart_ptr on multx_array array
    spv_float multy_array;	  //!< smart_ptr on multy_array array
    spv_float multz_array;	  //!< smart_ptr on multz_array array
  };
#endif //
