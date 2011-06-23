/*!
  \file ds_array_pool.cpp
  \brief Contains class for array

  Class #idata methods.

*/
#include "bs_bos_core_data_storage_stdafx.h"
#include "pool.h"

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
  : array_map(give_kernel::Instance().create_object(container::bs_type())) //ni(0), nj(0), nk(0),
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

  // create object
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (bs_array_map, (class)(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (bs_array_map, (class)(class));

  //BLUE_SKY_TYPE_IMPL_T (bos_array <base_strategy_di::item_array_t>, objbase, "bos_array_double", "", "");
  //BLUE_SKY_TYPE_IMPL_T (bos_array <base_strategy_fi::item_array_t>, objbase, "bos_array_float", "", "");
  //BLUE_SKY_TYPE_IMPL_T (bos_array <base_strategy_di::index_array_t>, objbase, "bos_array_int", "", "");

  typedef unsigned char uint8_t;
  typedef float float16_t;
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <strategy_t::index_t, uint8_t>),   1, (objbase), "bs_array_map_uint8",   "", "", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <strategy_t::index_t, float16_t>), 1, (objbase), "bs_array_map_float16", "", "", false);

  /*
  	\brief register plugin
   */
  bool register_array_map (const plugin_descriptor& pd)
  {
    return BS_KERNEL.register_type (pd, bs_array_map <strategy_t::index_t, uint8_t>::bs_type ())
      && BS_KERNEL.register_type (pd, bs_array_map <strategy_t::index_t, float16_t>::bs_type ())
      ;
  }

  template< > type_descriptor 
  bos_val_table <int, bs_array_map <strategy_t::index_t, uint8_t>::array_info >::bs_type()
  {
    return td_maker(std::string("_int_map_uint8"));
  }

  template< > type_descriptor 
  bos_val_table <int, bs_array_map <strategy_t::index_t, float16_t>::array_info >::bs_type()
  {
    return td_maker(std::string("_int_map_float16"));
  }


  //template< > type_descriptor bos_val_table< int, sp_obj >::bs_type()
  //{
  //  return td_maker(std::string("_int_sp_obj"));
  //}
};
