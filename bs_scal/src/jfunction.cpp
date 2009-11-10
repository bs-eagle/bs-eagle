/**
 * \file jfunction.cpp
 * \brief BLUE_SKY_TYPE_DECL impl
 * \author Sergey Miryanov
 * \date 22.05.2008
 * */
#include "bs_scal_stdafx.h"
#include "jfunction.h"

namespace blue_sky
  {

  template <typename strategy_t>
  jfunction<strategy_t>::jfunction(blue_sky::bs_type_ctor_param param)
      : is_valid (false)
      , perm_type (JFUNC_PERM_XY)
  {
    init (perm_type);
    is_valid = false;
  }
  template <typename strategy_t>
  jfunction<strategy_t>::jfunction(const this_t &j)
  : bs_refcounter (j), objbase (j)
  {
    st_phase    = j.st_phase;
    alpha       = j.alpha;
    beta        = j.beta;
    is_valid    = j.is_valid;
    plane_a     = j.plane_a;
    plane_b     = j.plane_b;
    perm_type   = j.perm_type;
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(jfunction,(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(jfunction,(class));
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jfunction<base_strategy_fi>), 1, (objbase), "jfunction_fi", "jfunction for capillary calculation fi", "jfunction for capillary calculation fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jfunction<base_strategy_di>), 1, (objbase), "jfunction_di", "jfunction for capillary calculation di", "jfunction for capillary calculation di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jfunction<base_strategy_mixi>), 1, (objbase), "jfunction_mixi", "jfunction for capillary calculation mixi", "jfunction for capillary calculation mixi", false);


} // namespace blue_sky
