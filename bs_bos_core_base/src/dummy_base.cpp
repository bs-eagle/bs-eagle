/**
 * \file dummy_base.cpp
 * \brief impl of
 * \author Sergey Miryanov
 * \date 30.03.2009
 * */
#include "bs_bos_core_base_stdafx.h"
#include "dummy_base.h"

namespace blue_sky {

  //////////////////////////////////////////////////////////////////////////
  // dummy linear solver base

  //! constructor
  dummy_base::dummy_base (bs_type_ctor_param /*param  = NULL */)
  {

  }

  //! copy constructor
  dummy_base::dummy_base (const dummy_base &d) : bs_refcounter (), objbase ()
  {
    if (this != &d)
      *this = d;
  }
  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE(dummy_base);
  BLUE_SKY_TYPE_STD_COPY(dummy_base);
  BLUE_SKY_TYPE_IMPL(dummy_base, objbase, "dummy_base", "Dummy base for python wrappers", "Dummy base for python wrappers!");


}
