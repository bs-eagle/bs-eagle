#ifndef FLUX_CONNECTIONS_IFACE_H
#define FLUX_CONNECTIONS_IFACE_H
/*! \file flux_connections_iface.h
	\brief This file declare interface class for bs flux connection data
	\author Mark Khait
	\date 2009-07-22
 * */
 
#ifndef PURE_MESH
      #include "../../bs_mtx/include/bcsr_matrix_iface.h"
#endif


namespace blue_sky
{
  
  class BS_API_PLUGIN flux_connections_iface : public objbase
  {
  ///////////////////////////////
  //  INTERNAL TYPE DECLARATION
  ///////////////////////////////
  public:
    
#ifndef PURE_MESH
      typedef smart_ptr <bcsr_matrix_iface, true>           sp_bcsr_t;
#else
      typedef csr_matrix *          sp_bcsr_t;
#endif

    //-------------------------------------------
    //  METHODS
    //===========================================
  public:

    //! default destructor
    virtual ~flux_connections_iface () {};
    
    //! get transmissibility matrix
    virtual sp_bcsr_t get_conn_trans () = 0;
    
    //! get matrix pointers in positive direction 
    virtual spv_long get_matrix_block_idx_plus () = 0;
    
    //! get matrix pointers in negative direction 
    virtual spv_long get_matrix_block_idx_minus () = 0;
  };
  
}; //namespace blue_sky
#endif //flux_connections_iface
