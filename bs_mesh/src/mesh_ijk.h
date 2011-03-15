#ifndef MESH_IJK_H
#define MESH_IJK_H
/**
	\file mesh_grdecl.h
 	\brief This file declare class for working with standart ijk mesh
	\author Iskhakov Ruslan
	\date 2008-11-13 */

#include "rs_smesh_base.h"
#include "mesh_element3d.h"
#include "flux_connections_iface.h"

//! class for basic work with mesh based on IJK and tpfa calculating

class  mesh_ijk : public rs_smesh_base
  {
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
    typedef grd_ecl::fpoint3d                           fpoint3d;
    typedef std::vector<grd_ecl::fpoint2d>              g_fpoint2d_vector;
    typedef boost::array <t_double, 3>                    center_t;

//-------------------------------------------
//  METHODS
//===========================================

  public:
    //! default constructor
    mesh_ijk();

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
    int build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn,
                                             spv_long boundary_array);

    /*!	\brief  find neighbours (adjacency matrix)
    	\param neig_matrix - bcsr adjacency matrix
    	\return 0 if success */
    int find_neighbours(sp_bcsr_t /*neig_matrix*/);

    /*! \brief get_block_dx_dy_dz
        \param n_elem - block number
        \param (return value) dx, dy, dz - size of current block*/
    void get_block_dx_dy_dz(t_long n_elem, t_double &dx, t_double &dy, t_double &dz)const;

    t_double get_block_dx(t_long n_elem) const
      {
        return dx_array[n_elem];
      };
    t_double get_block_dy(t_long n_elem) const
      {
        return dy_array[n_elem];
      };
    t_double get_block_dz(t_long n_elem) const
      {
        return dz_array[n_elem];
      };
    t_double get_depth(t_long n_elem) const
      {
        return (*depths)[n_elem];
      };
    t_double get_dtop(t_long n_elem) const
      {
        return (*depths)[n_elem];
      };

    //! get vertexes of cube (index = i+j*nx+k*nx*ny)
    //! length of (fpoint3d_vector) = 8
    
    grd_ecl::fpoint3d_vector calc_element (const t_long i, const t_long j, const t_long k) const;

    grd_ecl::fpoint3d_vector calc_element (const t_long index) const;
    
   
    //! get element center
    center_t
    get_center (t_long n_block) const;

    center_t
    get_center (t_long i, t_long j, t_long k) const;

  protected:
    /*! \brief calc depth for mesh (length n_active_elements)
        \return 0 if success*/
    int calc_depths();

  private:

    /*! \brief add index info into cols_ind, change m&p_memory block indexation, calculate and save transmissibility*/
    void inline set_neigbour_data (const t_long index1, const t_long index1_ext, const t_long index2_ext, t_long &conn_idx,
                        const t_long *rows_ptr, t_long *cols_ind, stdv_long &tmp_rows_ptr,
                        t_long *m_memory, t_long *p_memory,
                        t_long *cols_ind_tran, t_float *values_tran, direction dir);

    /*! \brief calculate transmissibility between two blocks
    	\param i - global indexation (X->Y->Z) for first block
    	\param j - global indexation (X->Y->Z) for second block (neigbor in positive direction)
    	\param d_dir - direction of transmissibility
    	\return transmissiblity value between 2 blocks */
    t_float calculate_tran(const t_long i, const t_long j, const direction d_dir) const;

    /*! \brief splicing block if it's too small
        \param volumes_temp - array of volumes
    	  \return number of splicing blocks  */
    int splicing(stdv_float & volumes_temp);

    //! calc values of all shift arrays
    int calc_shift_arrays();

//-------------------------------------------
//  VARIABLES
//===========================================
  protected:


    t_float *dx_array, *dy_array, *dz_array; //geometric properties
    t_float *tops_array; // layers

    stdv_float dx_shift_array, dy_shift_array, dz_shift_array; //array of shift for each directions
  }; //class mesh_ijk

#endif //MESH_IJK_H
