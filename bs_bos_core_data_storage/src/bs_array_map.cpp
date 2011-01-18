/*!
  \file ds_array_pool.cpp
  \brief Contains class for array

  Class #idata methods.

*/
#include "bs_bos_core_data_storage_stdafx.h"
#include "bs_array_map.h"

namespace blue_sky
  {
  //! empty destructor
  template <typename K, typename T>
  blue_sky::bs_array_map <K, T>::~bs_array_map () 
  {
  }


  //! default constructor
  template <typename K, typename T>
  bs_array_map <K, T>::bs_array_map (bs_type_ctor_param)
  : array_map(BS_KERNEL.create_object (container::bs_type())) //ni(0), nj(0), nk(0),
  {
  }

  //! copy constructor
  template<typename K, typename T>
  bs_array_map <K, T>::bs_array_map (const bs_array_map <K, T> &bs_array_map_ptr)
  : bs_refcounter (bs_array_map_ptr), objbase (bs_array_map_ptr),
  array_map(give_kernel::Instance().create_object_copy(bs_array_map_ptr.array_map))
  {
    ni = bs_array_map_ptr.ni;
    nj = bs_array_map_ptr.nj;
    nk = bs_array_map_ptr.nk;

    array_map->insert (bs_array_map_ptr.array_map->begin (), bs_array_map_ptr.array_map->end ());
  }

  /*
  	\brief initialization arrays parameters from table of values
    \param m_size - number of arrays
    \param dim - dimension elements
  */
  template <typename K, typename T>
  void  blue_sky::bs_array_map<K, T>::init (index_t nx, index_t ny, index_t nz)
  {
    if (nx == 0 || ny == 0 || nz == 0)
      throw bs_exception("bs_array_hash_map", "Invalid dimension sizes.");

    array_map->clear();

    ni = nx;
    nj = ny;
    nk = nz;
  }

  //bs stuff
  
  // create object
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (bs_array_map, (class)(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (bs_array_map, (class)(class));

  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <base_strategy_did::i_type_t, base_strategy_did::i_type_t>),   1, (objbase), "bs_array_map_ii",   "", "", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <base_strategy_dld::i_type_t, base_strategy_dld::i_type_t>),   1, (objbase), "bs_array_map_ll",   "", "", false);
  
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <base_strategy_did::i_type_t, base_strategy_did::fp_type_t>),   1, (objbase), "bs_array_map_id",   "", "", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <base_strategy_dld::i_type_t, base_strategy_dld::fp_type_t>),   1, (objbase), "bs_array_map_ld",   "", "", false);
  
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <base_strategy_fif::i_type_t, base_strategy_fif::fp_type_t>),   1, (objbase), "bs_array_map_if",   "", "", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <base_strategy_flf::i_type_t, base_strategy_flf::fp_type_t>),   1, (objbase), "bs_array_map_lf",   "", "", false);
  
  

  /*
  	\brief register plugin
   */
  bool register_array_map (const plugin_descriptor& pd)
  {
    return BS_KERNEL.register_type (pd, bs_array_map <base_strategy_did::i_type_t, base_strategy_did::i_type_t>::bs_type ())
      && BS_KERNEL.register_type (pd, bs_array_map <base_strategy_fif::i_type_t, base_strategy_fif::fp_type_t>::bs_type ())
      ;
  }
  
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <base_strategy_did::i_type_t, base_strategy_did::i_type_t>::array_info>), 1, (objbase), "bos_map<str, array_info_ii>", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <base_strategy_dld::i_type_t, base_strategy_dld::i_type_t>::array_info>), 1, (objbase), "bos_map<str, array_info_ll>", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <base_strategy_did::i_type_t, base_strategy_did::fp_type_t>::array_info>), 1, (objbase), "bos_map<str, array_info_id>", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <base_strategy_dld::i_type_t, base_strategy_dld::fp_type_t>::array_info>), 1, (objbase), "bos_map<str, array_info_ld>", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <base_strategy_fif::i_type_t, base_strategy_fif::fp_type_t>::array_info>), 1, (objbase), "bos_map<str, array_info_if>", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <base_strategy_flf::i_type_t, base_strategy_flf::fp_type_t>::array_info>), 1, (objbase), "bos_map<str, array_info_lf>", "map of str and facility", "map of str and facility", false);
  
  //BLUE_SKY_TYPE_IMPL_T_EXT (2, (bos_val_table <std::string, facility_manager<base_strategy_fi>::sp_facility_t>), 1, (objbase), "bos_map<str, facility>_seq_fi", "map of str and facility", "map of str and facility", false);
  /*
  template< > type_descriptor 
  bos_val_table <base_strategy_did::i_type_t, bs_array_map <base_strategy_did::i_type_t, base_strategy_did::i_type_t>::array_info >::bs_type()
  {
    return td_maker(std::string("bos_val_table_ii"));
  }
  
  template< > type_descriptor 
  bos_val_table <base_strategy_dld::i_type_t, bs_array_map <base_strategy_dld::i_type_t, base_strategy_dld::i_type_t>::array_info >::bs_type()
  {
    return td_maker(std::string("bos_val_table_ll"));
  }
  
  template< > type_descriptor 
  bos_val_table <base_strategy_did::i_type_t, bs_array_map <base_strategy_did::i_type_t, base_strategy_did::fp_type_t>::array_info >::bs_type()
  {
    return td_maker(std::string("bos_val_table_id"));
  }
  
  template< > type_descriptor 
  bos_val_table <base_strategy_dld::i_type_t, bs_array_map <base_strategy_dld::i_type_t, base_strategy_dld::fp_type_t>::array_info >::bs_type()
  {
    return td_maker(std::string("bos_val_table_ld"));
  }
  
   template< > type_descriptor 
  bos_val_table <base_strategy_did::i_type_t, bs_array_map <base_strategy_fif::i_type_t, base_strategy_fif::fp_type_t>::array_info >::bs_type()
  {
    return td_maker(std::string("bos_val_table_if"));
  }
  
  template< > type_descriptor 
  bos_val_table <base_strategy_flf::i_type_t, bs_array_map <base_strategy_flf::i_type_t, base_strategy_flf::fp_type_t>::array_info >::bs_type()
  {
    return td_maker(std::string("bos_val_table_lf"));
  }
*/

  //template< > type_descriptor bos_val_table< int, sp_obj >::bs_type()
  //{
  //  return td_maker(std::string("_int_sp_obj"));
  //}
};
