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


template<class strategy_t>
class BS_API_PLUGIN rs_smesh_base : public rs_mesh_base<strategy_t>
  {


  //-----------------------------------------
  // TYPES
  //-----------------------------------------

  public:
    ///////////////////////
    // BASE TYPES
    ///////////////////////
    typedef rs_mesh_base <strategy_t>                   base_t;

    typedef typename base_t::i_type_t                   i_type_t;
    typedef typename base_t::fp_type_t                  fp_type_t;
    typedef typename base_t::fp_storage_type_t          fp_storage_type_t;

    typedef typename base_t::index_array_t              index_array_t;
    typedef typename base_t::item_array_t               item_array_t;
    
    typedef typename base_t::sp_fp_array_t              sp_fp_array_t;
    typedef typename base_t::sp_i_array_t               sp_i_array_t;
    typedef typename base_t::sp_fp_storage_array_t      sp_fp_storage_array_t;
    typedef typename base_t::sp_bcsr_t                  sp_bcsr_t;
    typedef typename base_t::sp_idata_t                 sp_idata_t;
    typedef typename base_t::sp_flux_conn_iface_t       sp_flux_conn_iface_t;

    ///////////////////////
    // OWN TYPES
    ///////////////////////

    typedef rs_smesh_base <strategy_t>                   this_t;
    typedef boost::array <i_type_t, 3>                   index_point3d_t;
    typedef std::pair<i_type_t, i_type_t>                elem_index;
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
    void init_props (const sp_idata_t &idata);

    //! check mesh data
    void check_data() const;

    //! get min&max xyz coordinates
    void get_min_max_xyz (fp_type_t &amax_x, fp_type_t &amin_x, fp_type_t &amax_y, fp_type_t &amin_y, fp_type_t &amax_z, fp_type_t &amin_z) const
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


    /*	\brief count how many activ blocks in layers, forming according direction
    	\param dir - direction of layers dividing
    	\param elem_in_layers - (return value) vector of active block's number in layer[]
    	\return number of layers
    	implementation non-dependece from type of mesh (just from actnum_array)		 */
    int get_elems_n_in_layers(const direction d_dir, index_array_t &elem_in_layers) const;

    ////////indexation functions

    //! return index of block[i,j,k] - same as XYZ_to_inside
    i_type_t get_n_block (const i_type_t i, const i_type_t j, const i_type_t k) const
      {
        return XYZ_to_inside (i, j, k);
      }

    //! convert global [i,j,k] to local indexing
    i_type_t inline XYZ_to_inside(const i_type_t i, const i_type_t j, const i_type_t k) const
      {
        if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz)
          return -1;
        return (*base_t::base_t::ext_to_int)[i + j*nx + k*nx*ny];
      }

    //! convert global [i+j*nx+k*nx*ny] to local indexing
    i_type_t XYZ_to_inside(const i_type_t i) const
      {
        return (*base_t::base_t::ext_to_int)[i];
      }

    //! convert local index to global [i,j,k]
    void inside_to_XYZ(const i_type_t index, i_type_t &i1, i_type_t &j1, i_type_t &k1)const;

    //! convert local index to global [i+j*nx+k*nx*ny]
    i_type_t inside_to_XYZ(const i_type_t index) const
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
    virtual void get_block_dx_dy_dz(i_type_t /*n_elem*/, fp_type_t & /*dx*/, fp_type_t & /*dy*/, fp_type_t & /*dz*/) const
      {
        BS_ASSERT (false && "PURE CALL");
      };

    //! short  analog for previous function
    virtual float get_block_dx(i_type_t /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

    virtual float get_block_dy(i_type_t /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

    virtual float get_block_dz(i_type_t /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

    virtual float get_dtop(i_type_t /*n_elem*/) const
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      };

    // TODO: WTF??!!!
    virtual int get_coord_size() const
    {
      BS_ASSERT (false && "PURE CALL");
      return 0;
    }
    virtual int get_zcorn_size() const
    {
      BS_ASSERT (false && "PURE CALL");
      return 0;
    }

    virtual boost::array <fp_type_t, 3>
    get_center (i_type_t /*i*/, i_type_t /*j*/, i_type_t /*k*/) const
      {
        BS_ASSERT (false && "PURE CALL");

        static boost::array <fp_type_t, 3> dummy;
        return dummy;
      }

    virtual boost::array <fp_type_t, 3>
    get_center (i_type_t /*n_element*/) const
      {
        BS_ASSERT (false && "PURE CALL");

        static boost::array <fp_type_t, 3> dummy;
        return dummy;
      }

//-------------------------------------------
//  VARIABLES
//===========================================
  protected:
    i_type_t nx, ny, nz;	//!< grid dimens

    fp_type_t min_x, max_x; //!< min&max x coordinates
    fp_type_t min_y, max_y; //!< min&max y coordinates
    fp_type_t min_z, max_z; //!< min&max z coordinates

    fp_storage_type_t *permx_array;	  //!< smart_ptr on permx_array array
    fp_storage_type_t *permy_array;	  //!< smart_ptr on permy_array array
    fp_storage_type_t *permz_array;	  //!< smart_ptr on permz_array array

    fp_storage_type_t *multx_array;	  //!< smart_ptr on multx_array array
    fp_storage_type_t *multy_array;	  //!< smart_ptr on multy_array array
    fp_storage_type_t *multz_array;	  //!< smart_ptr on multz_array array
  };
#endif //
