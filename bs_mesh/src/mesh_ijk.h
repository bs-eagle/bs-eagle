#ifndef MESH_IJK_H
#define MESH_IJK_H
/**
	\file mesh_grdecl.h
 	\brief This file declare class for working with standart ijk mesh
	\author Iskhakov Ruslan
	\date 2008-11-13 */

#include "rs_smesh_base.h"
#include "fpoint3d.h"
#include "flux_connections_iface.h"

//! class for basic work with mesh based on IJK and tpfa calculating
template<class strategy_t>
class  mesh_ijk : public mesh_rs<strategy_t>
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

    ///////////////////////
    // OWN TYPES
    ///////////////////////
    typedef grd_ecl::fpoint3d                           fpoint3d;
    typedef std::vector<grd_ecl::fpoint2d>              g_fpoint2d_vector;
    typedef typename strategy_t::rhs_item_array_t       rhs_item_array_t;
    typedef boost::array <item_t, 3>                    center_t;
    
//-------------------------------------------
//  METHODS
//===========================================

  public:
    //! default constructor
    mesh_ijk() {};

    //! default destructor
    ~mesh_ijk() {};

    
    
    

    //-------------------------------------------
    //  INHERITED FUNCTIONS
    //===========================================
    
    using base_t::init_int_to_ext;
    using base_t::XYZ_to_inside;
    
    //! init mesh properties
    void init_props (const sp_idata_t &data);
    
    //! make indexation - create proxy and non-proxy array using info of minpv
    //! return number of non active blocks (less than mpv)
    int init_ext_to_int();
    
    //! check mesh data
    void check_data () const;

    /*!	\brief allocate jacobian and fill structure
    		\param n_block_size matrix block size
    		\param jacobian jacobian (like adjacency matrix)
    		\param flux_conn connection with calculated transmissibility
    		\return 0 if success*/
    int build_jacobian_and_flux_connections (const sp_bcsr_t &jacobian, const sp_flux_conn_iface_t &flux_conn, 
                                             index_array_t &boundary_array);

    /*!	\brief  find neighbours (adjacency matrix)
    	\param neig_matrix - bcsr adjacency matrix
    	\return 0 if success */
    int find_neighbours(sp_bcsr_t &/*neig_matrix*/);

    /*! \brief get_block_dx_dy_dz
        \param n_elem - block number
        \param (return value) dx, dy, dz - size of current block*/
    void get_block_dx_dy_dz(index_t n_elem, item_t &dx, item_t &dy, item_t &dz)const;
    float get_block_dx(index_t n_elem) const
      {
        return sp_dx[n_elem];
      };
    float get_block_dy(index_t n_elem) const
      {
        return sp_dy[n_elem];
      };
    float get_block_dz(index_t n_elem) const
      {
        return sp_dz[n_elem];
      };
    float get_depth(index_t n_elem) const
      {
        return depths[n_elem];
      };  
    float get_dtop(index_t n_elem) const
      {
        return depths[n_elem];
      };

    //! get vertexes of cube (index = i+j*nx+k*nx*ny)
    //! length of (fpoint3d_vector) = 8
    grd_ecl::fpoint3d_vector top_cube (const index_t i, const index_t j, const index_t k) const;

    grd_ecl::fpoint3d_vector top_cube (const index_t index) const;

    
    //! get element center
    center_t
    get_center (index_t n_block) const;
    
    center_t
    get_center (index_t i, index_t j, index_t k) const {return base_t::get_center (i,j,k);};

  protected:
    /*! \brief calc depth for mesh (length n_active_elements)
        \return 0 if success*/
    int calc_depths();

  private:

    /*! \brief add index info into cols_ind, change m&p_memory block indexation, calculate and save transmissibility*/
    void inline set_neigbour_data (const index_t index1, const index_t index1_ext, const index_t index2_ext, index_t &conn_idx, 
                        const index_array_t &rows_ptr, index_array_t &cols_ind, index_array_t& tmp_rows_ptr, 
                        index_array_t& m_memory, index_array_t& p_memory,
                        index_array_t &cols_ind_tran, rhs_item_array_t &values_tran, direction dir);

    /*! \brief calculate transmissibility between two blocks
    	\param i - global indexation (X->Y->Z) for first block
    	\param j - global indexation (X->Y->Z) for second block (neigbor in positive direction)
    	\param d_dir - direction of transmissibility
    	\return transmissiblity value between 2 blocks */
    item_t calculate_tran(const index_t i, const index_t j, const direction d_dir) const;

    /*! \brief splicing block if it's too small
        \param volumes_temp - array of volumes
    	  \return number of splicing blocks  */
    int splicing(item_array_t& volumes_temp);

    //! calc values of all shift arrays
    int calc_shift_arrays();

//-------------------------------------------
//  VARIABLES
//===========================================
  protected:
    using base_t::depths;
    using base_t::ext_to_int;
    using base_t::int_to_ext;
    using base_t::nx;
    using base_t::ny;
    using base_t::nz;
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
    

    blue_sky::array_float16_t sp_dx, sp_dy, sp_dz; //geometric properties
    blue_sky::array_float16_t sp_tops; // layers

    item_array_t dx_shift_array, dy_shift_array, dz_shift_array; //array of shift for each directions
  }; //class mesh_ijk

#endif //MESH_IJK_H
