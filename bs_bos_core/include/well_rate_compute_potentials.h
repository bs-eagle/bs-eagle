/**
 *       \file  well_rate_compute_potentials.h
 *      \brief  Compute perforation potentials for well derivs
 *              computation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_COMPUTE_POTENTIALS_H_
#define BS_WELLS_WELL_RATE_CONTROL_COMPUTE_POTENTIALS_H_

namespace blue_sky
{

  template <typename data_t, typename params_t>
  inline void
  compute_potentials (const data_t &data, params_t &params)
  {
    typedef typename params_t::item_t item_t;

    const item_t &bhp = params.perf_bhp;
    const item_t &po  = params.pressure[params.n_block];

    //const item_t &ro  = params.rho;
    //const item_t &g   = params.gravity;
    //const item_t &h   = params.diff_depth;

    params.Po = bhp /*+ (ro * g * h)*/ - po;
    params.Pw = params.Po;
    params.Pg = params.Po;

    if (params.is_w)
      params.Pw -= CAP_PRESSURE (data, params.phase_d, FI_PHASE_WATER);

    if (params.is_g)
      params.Pg -= CAP_PRESSURE (data, params.phase_d, FI_PHASE_GAS);

#ifdef _DEBUG
    //BOSOUT (section::wells, level::debug) << boost::format ("[%d]: bhp: %.10e, rho: %.10e, h: %.10e, pref: %.10e") % params.n_block % bhp % ro % h % po << bs_end;
    BOSOUT (section::wells, level::debug) << boost::format ("[%d]: bhp: %.10e, pref: %.10e") % params.n_block % bhp % po << bs_end;
    BOSOUT (section::wells, level::debug) << boost::format ("[%d]: pw: %.10e pg: %.10e po: %.10e") % params.n_block % params.Pw % params.Pg % params.Po << bs_end;
#endif
  }

  template <typename params_t>
  inline bool
  is_prod_potential (const params_t &params)
  {
    return params.Po < typename params_t::item_t (0.0);
  }

} // namespace blue_sky



#endif  // #ifndef BS_WELLS_WELL_RATE_CONTROL_COMPUTE_POTENTIALS_H_

