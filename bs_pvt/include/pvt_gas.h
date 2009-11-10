/**
 * \file pvt_gas.h
 * \brief
 * \author Miryanov Sergey
 * \date 29.04.2008
 */
#ifndef BS_PVT_GAS_H_
#define BS_PVT_GAS_H_

#include "pvt_base.h"

namespace blue_sky
  {

  /**
   * \brief pvt_gas
   */
  template <typename strategy_t>
  class pvt_gas : public pvt_base <strategy_t>
    {
    public:

      typedef strategy_t                        pvt_strategy_t;
      typedef pvt_base <strategy_t>             base_t;
      typedef typename base_t::item_t           item_t;
      typedef typename base_t::index_t          index_t;
      typedef typename base_t::index_array_t    index_array_t;
      typedef typename base_t::item_array_t     item_array_t;
      typedef typename base_t::input_vector_t   input_vector_t;
      typedef typename base_t::vector_t         vector_t;

      /**
       * \brief store values into data
       *
       * \param seq_vector
       */
      virtual void insert_vector (const input_vector_t &vec);

      /**
       * \brief generate interpolated data
       */
      void build (item_t atm_p, item_t min_p, item_t max_p, index_t n_intervals);

      /**
       * \brief calculate interpolated value
       *
       * \param[in] p
       * \param[out] inv_fvf
       * \param[out] d_inv_fvf
       * \param[out] inv_visc
       * \param[out] d_inv_visc
       * \param[out] inv_visc_fvf
       * \param[out] d_inv_visc_fvf
       */
      virtual void calc (const item_t p, item_t  *inv_fvf, item_t *d_inv_fvf,
                         item_t *inv_visc, item_t *d_inv_visc,
                         item_t *inv_visc_fvf, item_t *d_inv_visc_fvf) const;

      virtual void
      print () const;

    private:

      void check_gas ();
      void check_gas_internal ();

    private:

      vector_t main_pressure_;
      vector_t main_fvf_;
      vector_t main_visc_;

      vector_t pressure_;
      vector_t inv_fvf_;
      vector_t inv_visc_;
      vector_t inv_visc_fvf_;

    public:

      BLUE_SKY_TYPE_DECL_T (pvt_gas);
    };

} // namespace blue_sky


#endif  // #ifndef BS_PVT_GAS_H_
