#ifndef RS_SMESH_IFACE_H
#define RS_SMESH_IFACE_H
/*!
	\file rs_smesh_iface.h
  \brief This file declares interface for bs reservoir simulation structured meshes
  \author Mark Khait
  \date 2009-07-20
 */

#include "rs_mesh_iface.h"
#include "fpoint3d.h"

namespace blue_sky
  {

  template<class strategy_t>
  class BS_API_PLUGIN rs_smesh_iface : virtual public rs_mesh_iface <strategy_t>
    {
//+++++++++++++++++++++++++++++++++++++++++++
//  INTERNAL TYPE DECLARATION
//===========================================
    public:
      ///////////////////////
      // BASE TYPES
      ///////////////////////
      typedef rs_mesh_iface <strategy_t>                  base_t;

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
      //////////////////////

      typedef boost::array <index_t, 3>                   index_point3d_t;

    public:

      //! default destructor
      virtual ~rs_smesh_iface ()	{};

      //! get mesh dimensions
      virtual index_point3d_t get_dimens () = 0;

      //! return coords of block vertexes by IJK indexes
      virtual grd_ecl::fpoint3d_vector top_cube (const index_t i, const index_t j, const index_t k) const = 0;

      //! return coords of block vertexes by n_block index
      virtual grd_ecl::fpoint3d_vector top_cube (const index_t index) const = 0;

      //! return center point of an element
      virtual point3d_t get_element_center (const index_t index) const = 0;

      //! return center point of an element
      virtual point3d_t get_element_center (const index_t i, const index_t j, const index_t k) const = 0;

      //! return I, J and K structured mesh coordinates of an element by internal number
      virtual void get_element_int_to_ijk (const index_t n_element, index_t &i, index_t &j, index_t &k) const = 0;

      //! return internal number of an element by I, J and K structured mesh coordinates
      virtual index_t get_element_ijk_to_int (const index_t i, const index_t j, const index_t k) const = 0;
    };

};//namespace blue_sky
#endif // RS_SMESH_IFACE_H
