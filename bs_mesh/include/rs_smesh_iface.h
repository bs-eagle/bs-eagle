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

  
  class BS_API_PLUGIN rs_smesh_iface : public rs_mesh_iface 
    {
//+++++++++++++++++++++++++++++++++++++++++++
//  INTERNAL TYPE DECLARATION
//===========================================
    public:
      ///////////////////////
      // BASE TYPES
      ///////////////////////
      typedef rs_mesh_iface                   base_t;

      ///////////////////////
      // OWN TYPES
      //////////////////////

      typedef boost::array <t_long, 3>                   index_point3d_t;

    public:

      //! default destructor
      virtual ~rs_smesh_iface ()	{};

      //! get mesh dimensions
      virtual index_point3d_t get_dimens () = 0;

      //! return coords of block vertexes by IJK indexes
      virtual grd_ecl::fpoint3d_vector calc_element (const t_long i, const t_long j, const t_long k) const = 0;

      //! return coords of block vertexes by n_block index
      virtual grd_ecl::fpoint3d_vector calc_element (const t_long index) const = 0;
      

      //! return center point of an element
      virtual point3d_t get_element_center (const t_long index) const = 0;

      //! return center point of an element
      virtual point3d_t get_element_center (const t_long i, const t_long j, const t_long k) const = 0;

      //! return I, J and K structured mesh coordinates of an element by internal number
      virtual void get_element_int_to_ijk (const t_long n_element, t_long &i, t_long &j, t_long &k) const = 0;

      //! return internal number of an element by I, J and K structured mesh coordinates
      virtual t_long get_element_ijk_to_int (const t_long i, const t_long j, const t_long k) const = 0;

	  virtual spv_double
      get_element_sizes (const t_long n_element) = 0;
    };

};//namespace blue_sky
#endif // RS_SMESH_IFACE_H
