/**
 * @file fip_results_storage.cpp
 * @brief fip results storage
 * @author Sayfullin Ilshat
 * @date 2009-02-05
 */
#include "stdafx.h"
#include "fip_results_storage.h"

namespace blue_sky
  {


  /**
   * @brief constructor
   */
  fip_results::fip_results ()
  {
    d_params.resize (FIP_D_PARAM_COUNT);
  }

  /**
   * @brief destructor
   */
  fip_results::~fip_results ()
  {}

  /**
   * @brief constructor
   */
  fip_results_storage::fip_results_storage (bs_type_ctor_param /* = NULL */)
  {}
  fip_results_storage::fip_results_storage (const fip_results_storage& src)
      : bs_refcounter (src), objbase (src)
  {
    *this = src;
  }
  /**
   * @brief destructor
   */
  fip_results_storage::~fip_results_storage ()
  {}

  class fip_results;

  BLUE_SKY_TYPE_STD_CREATE (fip_results_storage)
  BLUE_SKY_TYPE_STD_COPY   (fip_results_storage)
  BLUE_SKY_TYPE_IMPL_SHORT (fip_results_storage, objbase, "fip_results_storage")
} //blue_sky

