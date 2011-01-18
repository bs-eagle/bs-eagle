#ifndef BS_DATA_MANAGER_H
#define BS_DATA_MANAGER_H
/*!
\file data_manager.h
\brief class data_manager
\author Morozov Andrey
*/

#include "data_class.h"
#include "keyword_manager.h"

namespace blue_sky {

  template <class strategy_t>
  class BS_API_PLUGIN data_manager: public objbase
    {
    public:
      //TYPES
      //typedef typename strategy_t::matrix_t       							matrix_t;
      //typedef typename strategy_t::item_array_t   							item_array_t;
      //typedef typename strategy_t::index_array_t  							index_array_t;
      typedef typename strategy_t::fp_storage_type_t         			fp_storage_type_t;
      typedef typename strategy_t::i_type_t        							  i_type_t;

      typedef idata<strategy_t>                                 idata_t;
      typedef smart_ptr <idata_t, true>										      sp_idata ;
      typedef smart_ptr <FRead, true>													  sp_reader_t;
      typedef keyword_manager<strategy_t>                       keyword_manager_t;
      typedef smart_ptr <keyword_manager_t, true>								sp_km_t;

      //METHODS
      ~data_manager();

    public:
    
      // initialize data manager
      void init();
      
      // read keyword file 
      void read_keyword_file(const std::string filename);
      
      // ACCESS
      
      // keyword manager accessor
      sp_idata get_data () {return data;};
      
      // keyword manager accessor
      sp_reader_t get_reader () {return reader;};
      
      // keyword manager accessor
      sp_km_t get_keyword_manager () {return km;};
      
      //CHECK
      void check_arrays_for_inactive_blocks () const;
      //void update_geometry() const;
      

    public:
      BLUE_SKY_TYPE_DECL_T(data_manager)

    public:
      sp_idata        data;												//!< data storage
      sp_reader_t     reader;                     //!< parser
      sp_km_t         km;                         //!< keyword manager
      locale_keeper   lkeeper;
    };

} //ns bs
#endif //BS_DATA_MANAGER_H
