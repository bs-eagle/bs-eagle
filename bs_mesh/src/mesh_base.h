#ifndef MESH_BASE_H
#define MESH_BASE_H
/*!
	\file mesh_base.h
	\brief This file declares base class for meshes
	\author Mark Khait
	\date 2009-07-16
 */

#ifndef PURE_MESH
  #include "conf.h"
  #include BS_FORCE_PLUGIN_IMPORT()
  #include "hdm_iface.h"
  #include "bcsr_matrix_iface.h"
  #include BS_STOP_PLUGIN_IMPORT()

  using namespace blue_sky;
#else
  #include "pure_mesh.h"
#endif

  
  class BS_API_PLUGIN mesh_base
    {

    //-----------------------------------------
    // TYPES
    //-----------------------------------------
    
    public:
    
      ///////////////////////
      // OWN TYPES
      ///////////////////////

#ifndef PURE_MESH
      typedef smart_ptr <hdm_iface, true>    sp_hdm_t;
      typedef smart_ptr <bcsr_matrix_iface, true>           sp_bcsr_t;
#else
      typedef csr_matrix *          sp_bcsr_t;
#endif

    //-----------------------------------------
    //  METHODS
    //-----------------------------------------

    public:
      
      ///////////////////////
      // INIT
      ///////////////////////

      //! default constructor
      mesh_base ();

      //! default destructor
      virtual ~mesh_base ()	{};

      //! init mesh
      virtual void init_props (const sp_hdm_t /*hdm*/) = 0;

      //! initialize ext_to_int indexation
      virtual int init_ext_to_int() = 0;
      
      //! initialize int_to_ext indexation
      virtual int init_int_to_ext() = 0;
      
      //! check mesh data
      void check_data () const;

      ///////////////////////
      // ACCESS
      ///////////////////////
      
      //! return number of active mesh elements
      t_long 
      get_n_active_elements ()const
      {
        return n_active_elements;
      }
      
      //! return number of mesh elements
      t_long 
      get_n_elements ()const
      {
        return n_elements;
      }
      
      //! get const ext_to_int
      t_long convert_ext_to_int (const t_long n_element) const
      {
#ifndef PURE_MESH
        return (*ext_to_int)[n_element];
#else
        return ext_to_int[n_element];
#endif
      }
      
      //! get const int_to_ext
      t_long get_element_int_to_ext (const t_long n_element) const
      {
#ifndef PURE_MESH
        return (*int_to_ext)[n_element];
#else
        return int_to_ext[n_element];
#endif
      }

      //! get const int_to_ext
      const spv_long get_int_to_ext() const
      {
        return int_to_ext;
      }
      
      //! get const ext_to_int
      const spv_long get_ext_to_int() const
      {
        return ext_to_int;
      }
      
      //! get mesh elements volumes
      const spv_float get_volumes () const
      {
        return volumes;
      }

      //! get connection_number
      t_long 
      get_n_connections() const
      {
        return n_connections;
      }

      ///////////////////////
      // WORK
      ///////////////////////

      //!  find neighbors and put it in neighbour matrix (adjacency matrix)
      virtual int find_neighbours(sp_bcsr_t /*neig_matrix*/) = 0;

    //-----------------------------------------
    //  VARIABLES
    //-----------------------------------------

    protected:

      t_long n_elements;	      //!< number of elements
      t_long n_active_elements;	//!< number of active elements
      t_long n_connections; //!< connection number
      
      //! indexations arrays
      spv_long ext_to_int;
      spv_long int_to_ext;
      
      spv_float volumes; //!< elements volumes
    };

#endif //
