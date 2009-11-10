/**
 * \file pvt_oil.h
 * \brief
 * \author Miryanov Sergey
 * \date 29.04.2008
 */
#ifndef BS_PVT_OIL_H_
#define BS_PVT_OIL_H_

#include "pvt_dead_oil.h"

namespace blue_sky
  {

  /**
   * \brief pvt_oil
   */
  template <typename strategy_t>
  class pvt_oil : public pvt_dead_oil <strategy_t>
    {
    public:

      typedef strategy_t                        pvt_strategy_t;
      typedef pvt_dead_oil <strategy_t>         base_t;
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
       *
       * \param amt_p atmospheric pressure
       * \param min_p minimal value of pressure
       * \param max_p maximal value of pressure
       * \param n_intervals number of intervals
       */
      virtual void build (item_t atm_p, item_t min_p, item_t max_p, int n_intervals);

      /**
       * \brief calculate interpolated value
       *
       * \param is_g
       * \param main_var
       * \param gor
       */
      virtual bool 
      calc (const bool is_g, const int main_var, const item_t p, const item_t gor,
            item_t *inv_fvf, item_t *d_inv_fvf, item_t *gor_d_inv_fvf,
            item_t *inv_visc, item_t *d_inv_visc, item_t *gor_d_inv_visc,
            item_t *inv_visc_fvf, item_t *d_inv_visc_fvf, item_t *gor_d_inv_visc_fvf,
            item_t *gas_oil_ratio, item_t *d_gas_oil_ratio,
            const item_t drsdt = -1.0, const item_t dt = 0,
            const item_t old_gas_oil_ratio = 0) const;

      virtual void
      print () const;

    private:

      /**
       * \brief check main data (PVTO in old terms) for physical reasons
       * \return
       */
      virtual void check_oil ();

      /**
       * \brief find respective interval of GOR
       *
       * \param gor gor for search
       * \param j1
       * \param j2
       * \param end_j1
       * \param end_j2
       */
      void get_compressibility_interval (item_t gor, index_t &j1, index_t &j2, index_t &end_j1, index_t &end_j2);

      /**
       * \brief generate interpolated data for compressibility
       *
       * \param n_intervals number of interpolation intervals
       */
      void build_compressibility (index_t n_intervals);

      /**
       * \brief calc undersaturated oil
       */
      inline bool 
      calc_undersaturated_oil (const item_t p, const item_t gor,
                               item_t *inv_fvf, item_t *d_inv_fvf, item_t *gor_d_inv_fvf,
                               item_t *inv_visc, item_t *d_inv_visc, item_t *gor_d_inv_visc,
                               item_t *inv_visc_fvf, item_t *d_inv_visc_fvf, item_t *gor_d_inv_visc_fvf,
                               item_t *gas_oil_ratio, item_t *d_gas_oil_ratio,
                               const item_t drsdt = -1.0, const item_t dt = 0,
                               const item_t old_gas_oil_ratio = 0) const;

    public:

      vector_t compress_fvf_;
      vector_t compress_visc_;

      using base_t::main_gpr_;
      using base_t::main_pressure_;
      using base_t::main_fvf_;
      using base_t::main_visc_;

      using base_t::pressure_;
      using base_t::gor_;
      using base_t::inv_fvf_;
      using base_t::inv_visc_;

    public:

      BLUE_SKY_TYPE_DECL_T (pvt_oil);

    };

} // namespace blue_sky

#endif  // #ifndef BS_PVT_OIL_H_
