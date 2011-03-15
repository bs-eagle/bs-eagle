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

  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <t_long, t_int>),   1, (objbase), "bs_array_map_ii",   "", "", false);
  //BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <t_long, t_long>),   1, (objbase), "bs_array_map_ll",   "", "", false);
  
  //BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <t_int, t_double>),   1, (objbase), "bs_array_map_id",   "", "", false);
  //BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <t_long, t_double>),   1, (objbase), "bs_array_map_ld",   "", "", false);
  
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <t_long, t_float>),   1, (objbase), "bs_array_map_if",   "", "", false);
  //BLUE_SKY_TYPE_IMPL_T_EXT (2, (bs_array_map <t_long, t_double>),   1, (objbase), "bs_array_map_lf",   "", "", false);
  
  

  /*
  	\brief register plugin
   */
  bool register_array_map (const plugin_descriptor& pd)
  {
    return BS_KERNEL.register_type (pd, bs_array_map <t_long, t_int>::bs_type ())
      && BS_KERNEL.register_type (pd, bs_array_map <t_long, t_float>::bs_type ())
      ;
  }
  
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <t_long, t_int>::array_info>), 1, (objbase), "bos_map<str, array_info_ii>", "map of str and facility", "map of str and facility", false);
  //BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <t_long, t_long>::array_info>), 1, (objbase), "bos_map<str, array_info_ll>", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <t_long, t_float>::array_info>), 1, (objbase), "bos_map<str, array_info_if>", "map of str and facility", "map of str and facility", false);
  //BLUE_SKY_TYPE_IMPL_T_EXT (3, (bos_val_table <std::string, bs_array_map <t_long, t_double>::array_info>), 1, (objbase), "bos_map<str, array_info_lf>", "map of str and facility", "map of str and facility", false);
};
