#ifndef RS_MESH_IFACE_H
#define RS_MESH_IFACE_H
/*!
	\file rs_mesh_iface.h
  \brief This file declares base class for bs reservoir simulation meshes
  \author Mark Khait
  \date 2009-07-20
 */

#include "flux_connections_iface.h"
#include "bs_array.h"

namespace blue_sky
  {

  class FRead;

  class /*BS_API_PLUGIN*/ rs_mesh_iface : public bs_node
    {
//+++++++++++++++++++++++++++++++++++++++++++
//  INTERNAL TYPE DECLARATION
//===========================================
    public:
      ///////////////////////
      // OWN TYPES
      ///////////////////////
      
      typedef rs_mesh_iface                  this_t;
      typedef smart_ptr <this_t, true>                    sp_this_t;
        
      //typedef std::vector<t_long>                       index_array_t;
      
      typedef flux_connections_iface                      flux_conn_iface_t;
      typedef smart_ptr <flux_conn_iface_t, true>         sp_flux_conn_iface_t;
      
      typedef bcsr_matrix_iface                           csr_matrix_t;
      typedef smart_ptr <csr_matrix_t, true>              sp_bcsr_t;

      
      typedef smart_ptr <hydrodynamic_model_iface, true>    sp_hdm_t;
      
      typedef boost::array <t_float, 3>                 point3d_t;
      typedef smart_ptr <FRead, true>										  sp_reader;
      
    public:

      //! default destructor
      virtual ~rs_mesh_iface ()	{};
      
      //! init mesh
      virtual void init_props (const sp_hdm_t hdm) = 0;

      //! initialize int_to_ext indexation
      virtual int init_int_to_ext() = 0;

      //! initialize int_to_ext indexation
      virtual int init_ext_to_int() = 0;
      
      //! check mesh data
      virtual void check_data () const = 0;
      
      //! return number of active mesh elements
      virtual t_long get_n_active_elements ()const = 0;
      
      //! return number of mesh elements
      virtual t_long get_n_elements ()const = 0;
      
      //! return number of mesh elements connections
      virtual t_long get_n_connections ()const = 0;
      
      //! //! return mesh dimensions range
      virtual void get_dimensions_range (t_double &dim1_max, t_double &dim1_min,
                                         t_double &dim2_max, t_double &dim2_min,
                                         t_double &dim3_max, t_double &dim3_min) const = 0;
      
      //! get const int_to_ext
      virtual const spv_long get_int_to_ext() const = 0;
      
      //! get const ext_to_int
      virtual const spv_long get_ext_to_int() const = 0;
      
      //! get mesh elements volumes
      virtual const spv_float get_volumes () const = 0;
      
      
      //! set darcy constant for correct transmissibility calculation 
      virtual void set_darcy (double darcy_constant_) = 0;
      
      //! get element internal number by external
      virtual t_long convert_ext_to_int (const t_long n_element) const = 0;

      //! get element external number by internal
      virtual t_long get_element_int_to_ext (const t_long n_element) const = 0;
      
      //! return element size in all 3 dimensions
      virtual void get_element_size (const t_long n_element, t_double &d_dim1, t_double &d_dim2, t_double &d_dim3) const = 0;

      //! return element size in 3rd dimension
      virtual t_double get_element_dim3_size (const t_long n_element) const = 0;
      
      //! return depth of mesh element
      virtual t_double get_element_depth(const t_long n_element) const = 0;
        
      //! return number of mesh elements connections
      virtual t_float get_element_dtop(const t_long n_element) const = 0;
      
      //! return center point of an element
      virtual point3d_t get_element_center (const t_long n_element)const = 0;
      
      //! return depths of cell centers (length n_active_elements)
      virtual const spv_float get_depths () const = 0;
      
      //! allocate jacobian 
      virtual int build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn, spv_long boundary_array) = 0;

    };

};//namespace blue_sky
#endif // RS_MESH_IFACE_H
