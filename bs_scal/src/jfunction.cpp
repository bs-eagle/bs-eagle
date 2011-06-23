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

  jfunction::jfunction(blue_sky::bs_type_ctor_param param)
      : is_valid (false)
      , perm_type (JFUNC_PERM_XY)
  {
    init (perm_type);
    is_valid = false;
  }
  jfunction::jfunction(const this_t &j)
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
  BLUE_SKY_TYPE_STD_CREATE (jfunction);
  BLUE_SKY_TYPE_STD_COPY (jfunction);
  BLUE_SKY_TYPE_IMPL (jfunction, objbase, "jfunction", "jfunction", "jfunction");


} // namespace blue_sky
