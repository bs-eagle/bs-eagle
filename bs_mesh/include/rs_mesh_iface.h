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

  template <class strategy_t> class idata;
  class FRead;

  template<class strategy_t>
  class /*BS_API_PLUGIN*/ rs_mesh_iface : public bs_node
    {
//+++++++++++++++++++++++++++++++++++++++++++
//  INTERNAL TYPE DECLARATION
//===========================================
    public:
      ///////////////////////
      // OWN TYPES
      ///////////////////////
      
      typedef strategy_t                                  strategy_type;
      
      typedef rs_mesh_iface <strategy_t>                  this_t;
      typedef smart_ptr <this_t, true>                    sp_this_t;
        
      typedef typename strategy_t::i_type_t               i_type_t;
      typedef typename strategy_t::fp_type_t              fp_type_t;
      
      //typedef std::vector<i_type_t>                       index_array_t;
      
      typedef smart_ptr<bs_array<fp_type_t>, true>        sp_fp_array_t;
      typedef smart_ptr<bs_array<i_type_t>, true>         sp_i_array_t;
            
      typedef flux_connections_iface<strategy_t>          flux_conn_iface_t;
      typedef smart_ptr <flux_conn_iface_t, true>         sp_flux_conn_iface_t;
      
      typedef bcsr_matrix_iface<strategy_t>               csr_matrix_t;
      typedef smart_ptr <csr_matrix_t, true>              sp_bcsr_t;

      typedef idata<strategy_t>                           idata_t;
      typedef smart_ptr <idata_t, true>                   sp_idata_t;
      
      typedef boost::array <fp_type_t, 3>                 point3d_t;
      typedef smart_ptr <FRead, true>										  sp_reader;
      
    public:

      //! default destructor
      virtual ~rs_mesh_iface ()	{};
      
      //! init mesh
      virtual void init_props (const sp_idata_t &idata) = 0;

      //! initialize int_to_ext indexation
      virtual int init_int_to_ext() = 0;

      //! initialize int_to_ext indexation
      virtual int init_ext_to_int() = 0;
      
      //! check mesh data
      virtual void check_data () const = 0;
      
      //! return number of active mesh elements
      virtual i_type_t get_n_active_elements ()const = 0;
      
      //! return number of mesh elements
      virtual i_type_t get_n_elements ()const = 0;
      
      //! return number of mesh elements connections
      virtual i_type_t get_n_connections ()const = 0;
      
      //! //! return mesh dimensions range
      virtual void get_dimensions_range (fp_type_t &dim1_max, fp_type_t &dim1_min,
                                         fp_type_t &dim2_max, fp_type_t &dim2_min,
                                         fp_type_t &dim3_max, fp_type_t &dim3_min) const = 0;
      
      //! get const int_to_ext
      virtual const sp_i_array_t get_int_to_ext() const = 0;
      
      //! get const ext_to_int
      virtual const sp_i_array_t get_ext_to_int() const = 0;
      
      //! get mesh elements volumes
      virtual const sp_fp_array_t get_volumes () const = 0;
      
      
      //! set darcy constant for correct transmissibility calculation 
      virtual void set_darcy (double darcy_constant_) = 0;
      
      //! get element internal number by external
      virtual i_type_t convert_ext_to_int (const i_type_t n_element) const = 0;

      //! get element external number by internal
      virtual i_type_t get_element_int_to_ext (const i_type_t n_element) const = 0;
      
      //! return element size in all 3 dimensions
      virtual void get_element_size (const i_type_t n_element, fp_type_t &d_dim1, fp_type_t &d_dim2, fp_type_t &d_dim3) const = 0;

      //! return element size in 3rd dimension
      virtual fp_type_t get_element_dim3_size (const i_type_t n_element) const = 0;
      
      //! return depth of mesh element
      virtual fp_type_t get_element_depth(const i_type_t n_element) const = 0;
        
      //! return number of mesh elements connections
      virtual fp_type_t get_element_dtop(const i_type_t n_element) const = 0;
      
      //! return center point of an element
      virtual point3d_t get_element_center (const i_type_t n_element)const = 0;
      
      //! return depths of cell centers (length n_active_elements)
      virtual const sp_fp_array_t get_depths () const = 0;
      
      //! allocate jacobian 
      virtual int build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn, sp_i_array_t boundary_array) = 0;

    };

};//namespace blue_sky
#endif // RS_MESH_IFACE_H
