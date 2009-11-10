/**
 * \file facility_manager.h
 * \brief facilities manager
 * \author Sergey Miryanov
 * \date 16.07.2008
 * */
#ifndef BS_FACILITY_MANAGER_H_
#define BS_FACILITY_MANAGER_H_

#include "facility_base.h"
#include "well_iterator.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class well;

  class data_storage_interface;

  template <typename strategy_t>
  class BS_API_PLUGIN facility_manager : public objbase
    {
    public:

      typedef facility_base <strategy_t>										facility_base_t;
      typedef well <strategy_t>															well_t;

      typedef smart_ptr <facility_base_t>										sp_facility_t;
      typedef bos_val_table <std::string, sp_facility_t>		facility_map_t;
      typedef smart_ptr <facility_map_t>										sp_facility_map_t;
      typedef smart_ptr <well_t, true>											sp_well_t;

      typedef smart_ptr <data_storage_interface, true>			sp_storage_t;

      typedef typename facility_map_t::iterator             facility_iterator_t;
      typedef typename facility_map_t::const_iterator       facility_const_iterator_t;

      //typedef well_iterator <locked_facility_map_t, facility_iterator_t, locked_well_t>  well_iterator_t;
      typedef facility_iterator_t                           well_iterator_t;
      typedef facility_const_iterator_t                     well_const_iterator_t;

    public:

      sp_well_t get_well (const std::string &group_name, const std::string &well_name) const;
      sp_well_t get_well (const std::string &well_name) const;

      void add_well (const sp_well_t &well);

      void save_data (const sp_storage_t &storage) const;

      well_const_iterator_t wells_begin () const;
      well_const_iterator_t wells_end () const;
      //well_iterator_t wells_begin ();
      //well_iterator_t wells_end ();

      //unsigned int wells_size () const;

      BLUE_SKY_TYPE_DECL_T (facility_manager);

    private:

      sp_facility_map_t					facility_list_;

    };

  bool
  facility_manager_register_type (const blue_sky::plugin_descriptor &pd);


}	// namespace blue_sky


#endif	// #ifndef BS_FACILITY_MANAGER_H_

