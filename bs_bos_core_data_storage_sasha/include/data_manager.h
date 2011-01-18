#ifndef BS_DATA_MANAGER_H
#define BS_DATA_MANAGER_H
/*!
\file data_manager.h
\brief class data_manager
\author Morozov Andrey
*/

#include "data_class.h"

namespace blue_sky {

  template <class strategy_t>
  class BS_API_PLUGIN data_manager: public objbase
    {
    public:
      //TYPES
      //typedef typename strategy_t::matrix_t       							matrix_t;
      //typedef typename strategy_t::item_array_t   							item_array_t;
      //typedef typename strategy_t::index_array_t  							index_array_t;
      typedef typename strategy_t::fp_type_t         							item_t;
      typedef typename strategy_t::i_type_t        							index_t;

      typedef idata <strategy_t>                                idata_t;
      typedef smart_ptr <idata_t, true>										      sp_idata ;
      typedef smart_ptr <FRead, true>													  sp_reader_t;

      //METHODS
      ~data_manager();

    public:
      //CHECK
      void check_arrays_for_inactive_blocks () const;
      //void update_geometry() const;

    public:
      BLUE_SKY_TYPE_DECL_T(data_manager)

    public:
      sp_idata        data;												//!< data storage
      sp_reader_t       reader;                     //!< parser
      locale_keeper   lkeeper;
    };

} //ns bs
#endif //BS_DATA_MANAGER_H
