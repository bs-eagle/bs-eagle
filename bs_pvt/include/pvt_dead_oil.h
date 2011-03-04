/**
 * \file pvt_dead_oil.h
 * \brief
 * \author Miryanov Sergey
 * \date 29.04.2008
 */
#ifndef BS_PVT_DEAD_OIL_H_
#define BS_PVT_DEAD_OIL_H_

#include "pvt_base.h"

namespace blue_sky
  {

  /**
   * \brief
   */
  class pvt_dead_oil : public pvt_base 
    {
    public:

      typedef pvt_base                          base_t;
      typedef t_double                          item_t;
      typedef t_int                             index_t;
      typedef v_double                          input_vector_t;
      typedef base_t::vector_t                  vector_t;

      /**
       * \brief parse line of chars and store values into data
       *
       * \param char_line line of chars
       * \return
       */
      //virtual void parse_char_line (const char * char_line);

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

      virtual
      item_t get_gor_for_pressure (item_t pressure_data) const;

      virtual
      item_t interpolate_and_fix (item_t cell_pbub) const;

      const vector_t &
      get_pressure () const
      {
        return pressure_;
      }
      const vector_t &
      get_gor () const
      {
        return gor_;
      }

      virtual void
      print () const;

    protected:

      /**
       * \brief check main data (PVDG(?) in old terms) for physical reasons
       * \return
       */
      virtual void check_oil ();

      /**
       * \brief generate interpolated data
       *
       * \param amt_p atmospheric pressure
       * \param min_p minimal value of pressure
       * \param max_p maximal value of pressure
       * \param n_intervals number of intervals
       * \param is_pvto is called for pvt_oil object or for pvt_dead_oil
       */
      int build_internal (item_t atm_p, item_t min_p, item_t max_p, int n_intervals, bool is_pvto);

      /**
       * \brief calc saturated oil
       */
      bool
      calc_saturated_oil (const bool is_g, const int main_var, const item_t p, const item_t gor,
                          item_t *inv_fvf, item_t *d_inv_fvf, item_t *gor_d_inv_fvf,
                          item_t *inv_visc, item_t *d_inv_visc, item_t *gor_d_inv_visc,
                          item_t *inv_visc_fvf, item_t *d_inv_visc_fvf, item_t *gor_d_inv_visc_fvf,
                          item_t *gas_oil_ratio, item_t *d_gas_oil_ratio,
                          const item_t drsdt = -1.0, const item_t dt = 0,
                          const item_t old_gas_oil_ratio = 0) const;

    protected:


      vector_t  main_gpr_;
      vector_t  main_pressure_;
      vector_t  main_fvf_;
      vector_t  main_visc_;

      vector_t  pressure_;
      vector_t  inv_fvf_;
      vector_t  inv_visc_;
      vector_t  inv_visc_fvf_;
      vector_t  gor_;

    public:

      BLUE_SKY_TYPE_DECL_T (pvt_dead_oil);

    };

} // namespace blue_sky

#endif  // #ifndef BS_PVT_DEAD_OIL_H_
