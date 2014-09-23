#ifndef FLUX_CONNECTIONS_H
#define FLUX_CONNECTIONS_H
/*! \file flux_connections.h
	\brief This file declare interface class which transfers flux connection data to the reservoir simulation process
	\author Iskhakov Ruslan
	\date 2008-05-20
 * */

#include "bs_common.h"
#include "bcsr_matrix_iface.h"
using namespace blue_sky;
 

class BS_API_PLUGIN flux_connections
{
  //+++++++++++++++++++++++++++++++++++++++++++
  //  INTERNAL TYPE DECLARATION
  //===========================================
public:
  typedef bcsr_matrix_iface                           csr_matrix_t;
  typedef smart_ptr <csr_matrix_t, true>              sp_bcsr_t;

  //-------------------------------------------
  //  METHODS
  //===========================================
public:
  //! default constructor
  flux_connections ();

  //! default destructor
  ~flux_connections () {};
  
  //! get transmissibility matrix
  sp_bcsr_t get_conn_trans () 
    {return conn_trans;};

  //! get matrix pointers in positive direction 
  spv_long  get_matrix_block_idx_plus ()
    {return matrix_block_idx_plus;};

  //! get matrix pointers in negative direction 
  spv_long  get_matrix_block_idx_minus ()
    {
      return matrix_block_idx_minus;
    };

  //-------------------------------------------
  //  VARIABLES
  //===========================================
public:

  //! CSR matrix for storing connections and tranmissibilities
  sp_bcsr_t conn_trans;
  spv_long matrix_block_idx_plus, matrix_block_idx_minus;//!< index in matrix for plus and minus blocks
  
};
#endif //flux_connections
