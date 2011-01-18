/** 
 * @file dens_matrix_tools.cpp
 * @brief 
 * @date 2009-11-24
 */
#include <time.h>
#include <stdlib.h>
#include <map>

#include "bs_mtx_stdafx.h"
#include "dens_matrix_tools.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  template <class strat_t>
  dens_matrix_tools<strat_t>::dens_matrix_tools (bs_type_ctor_param) 
        : dens_matrix_tools_iface<strat_t> ()
    {
    }
  template <class strat_t>
  dens_matrix_tools <strat_t>::dens_matrix_tools (const dens_matrix_tools <strat_t>& /*src*/) : bs_refcounter ()
     {
     }

template <class strat_t> int
dens_matrix_tools<strat_t>::random_init (sp_dens_t matrix, 
                                         const i_type_t ni, 
                                         const i_type_t nj,
                                         const i_type_t block_size,
                                         const fp_type_t rand_value_dispersion 
                                         ) const
{
  if (matrix->init (ni, nj, block_size))
    return -2;

  fp_storage_type_t *v = &(*(matrix->get_values ()))[0];

  srand ((unsigned)time( NULL ));
  
  for (i_type_t i = 0; i < ni; ++i)
    {
      for (i_type_t j = 0; j < nj; ++j)
        {
          v[i * nj + j] = (fp_storage_type_t)rand () / (fp_storage_type_t)RAND_MAX * (fp_storage_type_t)rand_value_dispersion;
        }
    }
  return 0;
}
      
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(dens_matrix_tools, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(dens_matrix_tools, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix_tools<base_strategy_fif>), 1,  (dens_matrix_tools_iface <base_strategy_fif>), "dens_matrix_tools_fif", "Tools for Dens Matrix class", "Tools realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix_tools<base_strategy_did>), 1,  (dens_matrix_tools_iface <base_strategy_did>), "dens_matrix_tools_did", "Tools for Dens Matrix class", "Tools realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix_tools<base_strategy_dif>), 1,  (dens_matrix_tools_iface <base_strategy_dif>), "dens_matrix_tools_dif", "Tools for Dens Matrix class", "Tools realization of Dens Matricies", false);

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix_tools<base_strategy_flf>), 1,  (dens_matrix_tools_iface <base_strategy_flf>), "dens_matrix_tools_flf", "Tools for Dens Matrix class", "Tools realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix_tools<base_strategy_dld>), 1,  (dens_matrix_tools_iface <base_strategy_dld>), "dens_matrix_tools_dld", "Tools for Dens Matrix class", "Tools realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix_tools<base_strategy_dlf>), 1,  (dens_matrix_tools_iface <base_strategy_dlf>), "dens_matrix_tools_dlf", "Tools for Dens Matrix class", "Tools realization of Dens Matricies", false);
}  // blue_sky namespace
