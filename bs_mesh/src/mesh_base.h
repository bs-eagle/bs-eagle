#ifndef MESH_BASE_H
#define MESH_BASE_H
/*!
	\file mesh_base.h
	\brief This file declares base class for meshes
	\author Mark Khait
	\date 2009-07-16
 */


using namespace blue_sky;

  template<class strategy_t>
  class BS_API_PLUGIN mesh_base
    {

    //-----------------------------------------
    // TYPES
    //-----------------------------------------
    
    public:
    
      ///////////////////////
      // OWN TYPES
      ///////////////////////
      typedef typename strategy_t::index_t                index_t;

      typedef typename strategy_t::index_array_t          index_array_t;
      typedef typename strategy_t::item_array_t           item_array_t;

      typedef idata                                       idata_t;
      typedef smart_ptr <idata_t, true>                   sp_idata_t;
      
      typedef typename strategy_t::csr_matrix_t           csr_matrix_t;
      typedef smart_ptr <csr_matrix_t, true>              sp_bcsr_t;

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
      virtual void init_props (const sp_idata_t &/*idata*/) = 0;
      
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
      index_t 
      get_n_active_elements ()const
      {
        return n_active_elements;
      }
      
      //! return number of mesh elements
      index_t 
      get_n_elements ()const
      {
        return n_elements;
      }
      
      //! get const ext_to_int
      index_t convert_ext_to_int (const index_t n_element) const
      {
        return ext_to_int[n_element];
      }
      
      //! get const int_to_ext
      index_t get_element_int_to_ext (const index_t n_element) const
      {
        return int_to_ext[n_element];
      }

      //! get const int_to_ext
      const index_array_t &
      get_int_to_ext() const
      {
        return int_to_ext;
      }
      
      //! get const ext_to_int
      const index_array_t &
      get_ext_to_int() const
      {
        return ext_to_int;
      }
      
      //! get mesh elements volumes
      const item_array_t &get_volumes () const
      {
        return volumes;
      }

      //! get connection_number
      index_t 
      get_n_connections() const
      {
        return n_connections;
      }

      ///////////////////////
      // WORK
      ///////////////////////

      //!  find neighbors and put it in neighbour matrix (adjacency matrix)
      virtual int find_neighbours(sp_bcsr_t &/*neig_matrix*/) = 0;

    //-----------------------------------------
    //  VARIABLES
    //-----------------------------------------

    protected:

      index_t n_elements;	      //!< number of elements
      index_t n_active_elements;	//!< number of active elements
      index_t n_connections; //!< connection number
      
      //! indexations arrays
      index_array_t ext_to_int;
      index_array_t int_to_ext;
      
      item_array_t volumes; //!< elements volumes
    };

#endif //
