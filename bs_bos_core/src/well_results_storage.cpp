/**
 * @file well_results_storage.cpp
 * @brief well results storage
 * @author Sayfullin Ilshat
 * @date 2009-01-28
 */
#include "stdafx.h"
#include "well_results_storage.h"

namespace blue_sky
  {

  /**
   * @brief constructor
   */
  connection_results::connection_results ()
  {
    d_params.resize (CONN_D_PARAM_COUNT);
    i_params.resize (CONN_I_PARAM_COUNT);
  }

  /**
   * @brief destructor
   */
  connection_results::~connection_results ()
  {}

  void
  connection_results::clear ()
  {
    for (size_t i = 0, cnt = d_params.size (); i < cnt; ++i)
      {
        d_params[i].clear ();
      }

    for (size_t i = 0, cnt = i_params.size (); i < cnt; ++i)
      {
        i_params[i].clear ();
      }
  }

  /**
   * @brief constructor
   */
  well_results::well_results ()
  {
    d_params.resize (WELL_D_PARAM_COUNT);
    i_params.resize (WELL_I_PARAM_COUNT);
  }

  /**
   * @brief destructor
   */
  well_results::~well_results ()
  {}

  void
  well_results::clear ()
  {
    dates.clear ();

    for (size_t i = 0, cnt = d_params.size (); i < cnt; ++i)
      {
        d_params[i].clear ();
      }

    for (size_t i = 0, cnt = i_params.size (); i < cnt; ++i)
      {
        i_params[i].clear ();
      }
  }

  /**
   * @brief constructor
   */
  well_results_storage::well_results_storage (bs_type_ctor_param param /* = NULL */)
  {}
  well_results_storage::well_results_storage (const well_results_storage& src)
      : bs_refcounter (src), objbase (src)
  {
    *this = src;
  }
  /**
   * @brief destructor
   */
  well_results_storage::~well_results_storage ()
  {}

  class well_results;
  class connection_results;

  BLUE_SKY_TYPE_STD_CREATE(well_results_storage)
  BLUE_SKY_TYPE_STD_COPY(well_results_storage)
  BLUE_SKY_TYPE_IMPL_SHORT(well_results_storage, objbase, "well_results_storage")
} //blue_sky

