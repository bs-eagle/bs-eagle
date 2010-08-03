/** 
 * @file jac_matrix.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-14
 */

#include "bs_matrix_stdafx.h"
#include "jac_matrix.h"

namespace blue_sky
  {

  //! constructor
  template <class strat_t>
  jac_matrix<strat_t>::jac_matrix (bs_type_ctor_param /*param*/)
      : sp_accum_matrix (BS_KERNEL.create_object (bdiag_matrix_t::bs_type ()))
      , sp_flux_matrix (BS_KERNEL.create_object (bcsr_matrix_t::bs_type ()))
      , sp_facility_matrix (BS_KERNEL.create_object (bcsr_matrix_t::bs_type ()))
  {

  }

  //! copy constructor
  template <class strat_t>
  jac_matrix<strat_t>::jac_matrix (const jac_matrix &matrix) 
        : bs_refcounter () 
  {
    BS_ASSERT (false && "TEST ME");

    if (this != &matrix)
      {
        *sp_accum_matrix = *(matrix.sp_accum_matrix);
        *sp_facility_matrix = *(matrix.sp_facility_matrix);
        *sp_flux_matrix = *(matrix.sp_flux_matrix);
      }

  }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(jac_matrix, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(jac_matrix, (class));


  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jac_matrix<base_strategy_fif>), 1,  (jac_matrix_iface <base_strategy_fif> ), "jac_matrix_fif", "Jacobian Matrix class", "Realization of Jacobian Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jac_matrix<base_strategy_did>), 1,  (jac_matrix_iface <base_strategy_did> ), "jac_matrix_did", "Jacobian Matrix class", "Realization of Jacobian Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jac_matrix<base_strategy_dif>), 1,  (jac_matrix_iface <base_strategy_dif> ), "jac_matrix_dif", "Jacobian Matrix class", "Realization of Jacobian Matricies", false);

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jac_matrix<base_strategy_flf>), 1,  (jac_matrix_iface <base_strategy_flf> ), "jac_matrix_flf", "Jacobian Matrix class", "Realization of Jacobian Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jac_matrix<base_strategy_dld>), 1,  (jac_matrix_iface <base_strategy_dld> ), "jac_matrix_dld", "Jacobian Matrix class", "Realization of Jacobian Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jac_matrix<base_strategy_dlf>), 1,  (jac_matrix_iface <base_strategy_dlf> ), "jac_matrix_dlf", "Jacobian Matrix class", "Realization of Jacobian Matricies", false);
  }


