#ifndef FLUX_CONNECTIONS_IFACE_H
#define FLUX_CONNECTIONS_IFACE_H
/*! \file flux_connections_iface.h
	\brief This file declare interface class for bs flux connection data
	\author Mark Khait
	\date 2009-07-22
 * */
 
#include "../../bs_mtx/include/bcsr_matrix_iface.h"

namespace blue_sky
{
  
   
  template<typename strategy_t>
  class BS_API_PLUGIN flux_connections_iface : virtual public objbase
  {
  ///////////////////////////////
  //  INTERNAL TYPE DECLARATION
  ///////////////////////////////
  public:
    typedef typename strategy_t::i_type_t               i_type_t;
    typedef smart_ptr <bs_array<i_type_t>, true>        sp_i_array_t;
    
    typedef bcsr_matrix_iface<strategy_t>               csr_matrix_t;
    typedef smart_ptr <csr_matrix_t, true>              sp_bcsr_t;

    //-------------------------------------------
    //  METHODS
    //===========================================
  public:

    //! default destructor
    virtual ~flux_connections_iface () {};
    
    //! get transmissibility matrix
    virtual sp_bcsr_t get_conn_trans () = 0;
    
    //! get matrix pointers in positive direction 
    virtual sp_i_array_t get_matrix_block_idx_plus () = 0;
    
    //! get matrix pointers in negative direction 
    virtual sp_i_array_t get_matrix_block_idx_minus () = 0;
  };
  
}; //namespace blue_sky
#endif //flux_connections_iface
