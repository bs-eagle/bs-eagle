#ifndef FLUX_CONNECTIONS_H
#define FLUX_CONNECTIONS_H
/*! \file flux_connections.h
	\brief This file declare interface class which transfers flux connection data to the reservoir simulation process
	\author Iskhakov Ruslan
	\date 2008-05-20
 * */

#include "strategies.h"
 
class flux_connections
{
  //+++++++++++++++++++++++++++++++++++++++++++
  //  INTERNAL TYPE DECLARATION
  //===========================================
public:
  typedef blue_sky::strategy_t                                  strategy_type;
  typedef blue_sky::strategy_t::index_array_t          index_array_t;

  typedef blue_sky::strategy_t::csr_matrix_t           csr_matrix_t;
  typedef blue_sky::smart_ptr <csr_matrix_t, true>    sp_bcsr_t;

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
  index_array_t &
  get_matrix_block_idx_plus ()
    {return matrix_block_idx_plus;};

  //! get matrix pointers in negative direction 
  index_array_t &
  get_matrix_block_idx_minus ()
    {
      return matrix_block_idx_minus;
    };

  //-------------------------------------------
  //  VARIABLES
  //===========================================
public:

  //! CSR matrix for storing connections and tranmissibilities
  sp_bcsr_t conn_trans;
  index_array_t matrix_block_idx_plus, matrix_block_idx_minus;//!< index in matrix for plus and minus blocks
  
};
#endif //flux_connections
