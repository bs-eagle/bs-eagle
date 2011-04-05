#ifndef BS_FLUX_CONNECTIONS_H
#define BS_FLUX_CONNECTIONS_H
/*! \file bs_flux_connections.h
	\brief This file declare bs wrapper over flux_connections
	\author Mark Khait
	\date 2009-07-22
 * */
 
#include "flux_connections.h"
#include "flux_connections_iface.h"
 
namespace blue_sky
  {
    
    class BS_API_PLUGIN bs_flux_connections : public flux_connections_iface
    {
      //+++++++++++++++++++++++++++++++++++++++++++
      //  INTERNAL TYPE DECLARATION
      //===========================================
    public:
      typedef flux_connections_iface  base_t;
      
      //-------------------------------------------
      //  METHODS
      //===========================================
    public:
      //! blue-sky class declaration
      BLUE_SKY_TYPE_DECL(bs_flux_connections);

      //! default destructor
      ~bs_flux_connections () {};
      
      //! get transmissibility matrix
      sp_bcsr_t get_conn_trans ()
        { return wrapped.get_conn_trans();};
      
      //! get matrix pointers in positive direction 
      spv_long get_matrix_block_idx_plus ()
        { return wrapped.get_matrix_block_idx_plus();};
      
      //! get matrix pointers in negative direction 
      spv_long get_matrix_block_idx_minus ()
        { return wrapped.get_matrix_block_idx_minus();};

    ////////////////////
    // wrapped class
    ///////////////////

    private:
      flux_connections wrapped;
      
    };
    
  }; //namespace blue_sky
#endif //bs_flux_connections
