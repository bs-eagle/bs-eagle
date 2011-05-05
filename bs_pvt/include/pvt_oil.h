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
  class BS_API_PLUGIN pvt_oil : public pvt_dead_oil 
    {
    public:

      typedef pvt_dead_oil                      base_t;
      typedef base_t::vector_t                  vector_t;

      enum {
         PVT_OIL_COMPRESS_FVF = PVT_OIL_TOTAL,
         PVT_OIL_COMPRESS_VISC,
         PVT_OIL_TOTAL_2
      };

      /**
       * \brief store values into data
       *
       * \param seq_vector
       */
      virtual void insert_vector (const v_double &vec);

      /**
       * \brief generate interpolated data
       *
       * \param amt_p atmospheric pressure
       * \param min_p minimal value of pressure
       * \param max_p maximal value of pressure
       * \param n_intervals number of intervals
       */
      virtual void build (t_double atm_p, t_double min_p, t_double max_p, t_long n_intervals);

      /**
       * \brief calculate interpolated value
       *
       * \param is_g
       * \param main_var
       * \param gor
       */
      virtual bool 
      calc (const bool is_g, const int main_var, const t_double p, const t_double gor,
            t_double *inv_fvf, t_double *d_inv_fvf, t_double *gor_d_inv_fvf,
            t_double *inv_visc, t_double *d_inv_visc, t_double *gor_d_inv_visc,
            t_double *inv_visc_fvf, t_double *d_inv_visc_fvf, t_double *gor_d_inv_visc_fvf,
            t_double *gas_oil_ratio, t_double *d_gas_oil_ratio,
            const t_double drsdt = -1.0, const t_double dt = 0,
            const t_double old_gas_oil_ratio = 0) const;

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
      void get_compressibility_interval (t_double gor, t_long &j1, t_long &j2, t_long &end_j1, t_long &end_j2);

      /**
       * \brief generate interpolated data for compressibility
       *
       * \param n_intervals number of interpolation intervals
       */
      void build_compressibility (t_long n_intervals);

      /**
       * \brief calc undersaturated oil
       */
      inline bool 
      calc_undersaturated_oil (const t_double p, const t_double gor,
                               t_double *inv_fvf, t_double *d_inv_fvf, t_double *gor_d_inv_fvf,
                               t_double *inv_visc, t_double *d_inv_visc, t_double *gor_d_inv_visc,
                               t_double *inv_visc_fvf, t_double *d_inv_visc_fvf, t_double *gor_d_inv_visc_fvf,
                               t_double *gas_oil_ratio, t_double *d_gas_oil_ratio,
                               const t_double drsdt = -1.0, const t_double dt = 0,
                               const t_double old_gas_oil_ratio = 0) const;

    public:

      //vector_t compress_fvf_;
      //vector_t compress_visc_;

      //using base_t::main_gpr_;
      //using base_t::main_pressure_;
      //using base_t::main_fvf_;
      //using base_t::main_visc_;

      //using base_t::pressure_;
      //using base_t::gor_;
      //using base_t::inv_fvf_;
      //using base_t::inv_visc_;

    public:

      BLUE_SKY_TYPE_DECL_T (pvt_oil);

    };

} // namespace blue_sky

#endif  // #ifndef BS_PVT_OIL_H_
