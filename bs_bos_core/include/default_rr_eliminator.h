/**
 *       \file  default_rr_eliminator.h
 *      \brief  Eliminates rr for default well implementation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  20.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Should be moved to src/ directory
 * */
#ifndef BS_DEFAULT_RR_ELIMINATOR_H_
#define BS_DEFAULT_RR_ELIMINATOR_H_

namespace blue_sky {

  /**
   * \class default_rr_eliminator
   * \brief Declaration of default_rr_eliminator
   * \todo  Draw formulas
   * */
  template <size_t size>
  struct default_rr_eliminator
  {
    /**
     * \brief  Eliminates diagonal element for well controlled by rate,
     *         res += (rr - rw * wr) * dt
     * */
    void 
    process_diag_rate ()
    {
    }
    /**
     * \brief  Eliminates diagonal element for well controlled by BHP,
     *         res += rr * dt
     * */
    void
    process_diag_bhp ()
    {
    }
    /**
     * \brief  Eliminates non-diagonal element,
     *         res += -(rw * wr) * dt
     * */
    void 
    process_rate ()
    {
    }
  };
  /**
   * \class default_rr_eliminator <3>
   * \brief Implements default_rr_eliminator for 3phase models
   * */
  template <>
  struct default_rr_eliminator <3>
  {
    template <typename res_t, typename rr_t, typename rw_t, typename wr_t, typename dt_t>
    static void
    process_diag_rate (res_t &res, const rr_t &rr, const rw_t &rw, const wr_t &wr, dt_t dt)
    {
      res[0] += (rr[0] - (rw[0] * wr[0])) * dt;
      res[1] += (rr[1] - (rw[0] * wr[1])) * dt;
      res[2] += (rr[2] - (rw[0] * wr[2])) * dt;
      res[3] += (rr[3] - (rw[1] * wr[0])) * dt;
      res[4] += (rr[4] - (rw[1] * wr[1])) * dt;
      res[5] += (rr[5] - (rw[1] * wr[2])) * dt;
      res[6] += (rr[6] - (rw[2] * wr[0])) * dt;
      res[7] += (rr[7] - (rw[2] * wr[1])) * dt;
      res[8] += (rr[8] - (rw[2] * wr[2])) * dt;
    }
    template <typename res_t, typename rr_t, typename dt_t>
    static void
    process_diag_bhp (res_t &res, const rr_t &rr, dt_t dt)
    {
      res[0] += rr[0] * dt;
      res[1] += rr[1] * dt;
      res[2] += rr[2] * dt;
      res[3] += rr[3] * dt;
      res[4] += rr[4] * dt;
      res[5] += rr[5] * dt;
      res[6] += rr[6] * dt;
      res[7] += rr[7] * dt;
      res[8] += rr[8] * dt;
    }
    template <typename res_t, typename rw_t, typename wr_t, typename dt_t>
    static void
    process_rate (res_t &res, const rw_t &rw, const wr_t &wr, dt_t dt)
    {
      res[0] += -(rw[0] * wr[0]) * dt;
      res[1] += -(rw[0] * wr[1]) * dt;
      res[2] += -(rw[0] * wr[2]) * dt;
      res[3] += -(rw[1] * wr[0]) * dt;
      res[4] += -(rw[1] * wr[1]) * dt;
      res[5] += -(rw[1] * wr[2]) * dt;
      res[6] += -(rw[2] * wr[0]) * dt;
      res[7] += -(rw[2] * wr[1]) * dt;
      res[8] += -(rw[2] * wr[2]) * dt;
    }
  };
  /**
   * \class default_rr_eliminator <2>
   * \brief Implements default_rr_eliminator for 2phase models
   * */
  template <>
  struct default_rr_eliminator <2>
  {
    template <typename res_t, typename rr_t, typename rw_t, typename wr_t, typename dt_t>
    static void
    process_diag_rate (res_t &res, const rr_t &rr, const rw_t &rw, const wr_t &wr, dt_t dt)
    {
      res[0] += (rr[0] - (rw[0] * wr[0])) * dt;
      res[1] += (rr[1] - (rw[0] * wr[1])) * dt;
      res[2] += (rr[2] - (rw[1] * wr[0])) * dt;
      res[3] += (rr[3] - (rw[1] * wr[1])) * dt;
    }
    template <typename res_t, typename rr_t, typename dt_t>
    static void
    process_diag_bhp (res_t &res, const rr_t &rr, dt_t dt)
    {
      res[0] += rr[0] * dt;
      res[1] += rr[1] * dt;
      res[2] += rr[2] * dt;
      res[3] += rr[3] * dt;
    }
    template <typename res_t, typename rw_t, typename wr_t, typename dt_t>
    static void
    process_rate (res_t &res, const rw_t &rw, const wr_t &wr, dt_t dt)
    {
      res[0] += -(rw[0] * wr[0]) * dt;
      res[1] += -(rw[0] * wr[1]) * dt;
      res[2] += -(rw[1] * wr[0]) * dt;
      res[3] += -(rw[1] * wr[1]) * dt;
    }
  };
  /**
   * \class default_rr_eliminator <1>
   * \brief Implements default_rr_eliminator for 1phase models
   * */
  template <>
  struct default_rr_eliminator <1>
  {
    template <typename res_t, typename rr_t, typename rw_t, typename wr_t, typename dt_t>
    static void
    process_diag_rate (res_t &res, const rr_t &rr, const rw_t &rw, const wr_t &wr, dt_t dt)
    {
      res[0] += (rr[0] - (rw[0] * wr[0])) * dt;
    }
    template <typename res_t, typename rr_t, typename dt_t>
    static void
    process_diag_bhp (res_t &res, const rr_t &rr, dt_t dt)
    {
      res[0] += rr[0] * dt;
    }
    template <typename res_t, typename rw_t, typename wr_t, typename dt_t>
    static void
    process_rate (res_t &res, const rw_t &rw, const wr_t &wr, dt_t dt)
    {
      res[0] += -(rw[0] * wr[0]) * dt;
    }
  };
} // namespace blue_sky

#endif // #ifndef BS_DEFAULT_RR_ELIMINATOR_H_

