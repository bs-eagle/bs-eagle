/**
 *       \file  norm_calc.cpp
 *      \brief  Calculates norms
 *     \author  Borschuk Oleg
 *       \date  31.01.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "norm_calc.h"

namespace blue_sky
  {

#ifdef _MPI
//#define _MPI_DEBUG_NORM
#endif //_MPI

  using namespace norms;

#define SET_VALUE(name_idx, value, descr, flag) \
  val[name_idx] = value;                        \
  name[name_idx] = descr;                       \
  p_flag[name_idx] = flag;

//! default constructor
  norms_storage::norms_storage ()
  {
    val.resize (NORMS_COUNTER);
    name.resize (NORMS_COUNTER);
    clear ();

    SET_VALUE (C_CPV, 0.0, "C-norm", 1);
    SET_VALUE (C_CPV_GAS, 0.0, "C-norm (gas)", 1);
    SET_VALUE (C_CPV_WATER, 0.0, "C-norm (water)", 1);
    SET_VALUE (C_CPV_OIL, 0.0, "C-norm (oil)", 1);

    SET_VALUE (C_ACPV, 0.0, "C-norm average", 0);
    SET_VALUE (C_ACPV_GAS, 0.0, "C-norm average (gas)", 0);
    SET_VALUE (C_ACPV_WATER, 0.0, "C-norm average (water)", 1);
    SET_VALUE (C_ACPV_OIL, 0.0, "C-norm average (oil)", 1);

    SET_VALUE (L2_CPV, 0.0, "L2-norm", 0);
    SET_VALUE (L2_CPV_GAS, 0.0, "L2-norm (gas)", 0);
    SET_VALUE (L2_CPV_WATER, 0.0, "L2-norm (water)", 0);
    SET_VALUE (L2_CPV_OIL, 0.0, "L2-norm (oil)", 0);

    SET_VALUE (L2_ACPV, 0.0, "L2-norm average", 0);
    SET_VALUE (L2_ACPV_GAS, 0.0, "L2-norm average (gas)", 0);
    SET_VALUE (L2_ACPV_WATER, 0.0, "L2-norm average (water)", 0);
    SET_VALUE (L2_ACPV_OIL, 0.0, "L2-norm average (oil)", 0);

    SET_VALUE (MB_ERR, 0.0, "Mat. Balance error", 1);
    SET_VALUE (MB_ERR_GAS, 0.0, "Mat. Balance error (gas)", 0);
    SET_VALUE (MB_ERR_WATER, 0.0, "Mat. Balance error (water)", 0);
    SET_VALUE (MB_ERR_OIL, 0.0, "Mat. Balance error (oil)", 0);

#ifdef _MPI_DEBUG_NORM //all norms will be printed
    memset (p_flag, -1, NORMS_COUNTER * sizeof (int));
#endif //_MPI
  }

//! default constructor
  norms_storage::~norms_storage ()
  {
  }

//! clear all norms
  void
  norms_storage::clear ()
  {
    val.assign (val.size (), 0);
    idx.assign (0);
  }

//! print norms
  /*void
  norms_storage::print (const fi_mesh *msh)
  {
    int n = val.size ();
    unsigned int i, j, k;

    rep->print (LOG_CHECK_SECTION, LOG_ERR, "========== NORM TYPE ==========|= RESIDUAL =|= CELL ==|=== I, J, K ===|\n", NORMS_COUNTER);
    for (int l = 0; l < n; ++l)
      {
        if (!p_flag[l])
          continue;
          msh->convert_num_to_block_coord (idx[l], i, j, k);
  #ifdef _MPI_DEBUG_NORM
          rep->switch_to_collective_print();
          rep->fflush_and_wait();
  #endif
          rep->print (LOG_CHECK_SECTION, LOG_ERR, "%30s |%11.2le |%8d |%4d,%4d,%4d |\n",
                    val[l].second.c_str (), val[l].first, idx[l], i + 1, j + 1, k + 1);
  #ifdef _MPI_DEBUG_NORM
          rep->switch_to_root_print();
          rep->fflush_and_wait();
  #endif
      }
    rep->print (LOG_CHECK_SECTION, LOG_ERR, "===============================|============|=========|===============|\n");
  }*/


  norms_storage &
  norms_storage::operator= (const norms_storage &rhs)
  {
    val = rhs.val;
    std::copy (rhs.idx.begin(),rhs.idx.end(), idx.begin());
    return *this;
  }

} //ns bs
