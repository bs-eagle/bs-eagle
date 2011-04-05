/**
 * \file pvt_water.h
 * \brief
 * \author Miryanov Sergey
 * \date 29.04.2008
 */
#ifndef BS_PVT_WATER_H_
#define BS_PVT_WATER_H_

#include "pvt_base.h"

namespace blue_sky
  {

  /**
   * \brief pvt_water
   */
  class pvt_water : public pvt_base 
    {
    public:

      typedef pvt_base                          base_t;
      typedef base_t::vector_t                  vector_t;
 
       enum {
         PVT_WATER_INPUT_PRESSURE = 0,
         PVT_WATER_INPUT_FVF,
         PVT_WATER_INPUT_COMPRESSIBILITY,
         PVT_WATER_INPUT_VISC,
         PVT_WATER_INPUT_VISCOSIBILITY,
         PVT_WATER_INPUT_GPR,
         PVT_WATER_INPUT_TOTAL
       };  

 
      enum {
         PVT_WATER_PRESSURE = 0,
         PVT_WATER_INV_FVF,
         PVT_WATER_INV_VISC,
         PVT_WATER_INV_VISC_FVF,
         PVT_WATER_TOTAL
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
      void build (t_double atm_p, t_double min_p, t_double max_p, t_long n_intervals);

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
      virtual void calc (const t_double p, t_double  *inv_fvf, t_double *d_inv_fvf,
                         t_double *inv_visc, t_double *d_inv_visc,
                         t_double *inv_visc_fvf, t_double *d_inv_visc_fvf) const;

      virtual void
      print () const;

    private:

      void check_water ();

    private:

      //vector_t main_pressure_;
      //vector_t main_compressibility_;
      //vector_t main_fvf_;
      //vector_t main_viscosibility_;
      //vector_t main_visc_;
      //vector_t main_gpr_;

      //vector_t pressure_;
      //vector_t inv_fvf_;
      //vector_t inv_visc_;
      //vector_t inv_visc_fvf_;

    public:

      BLUE_SKY_TYPE_DECL_T (pvt_water);

    };

} // namespace blue_sky


#endif  // #ifndef BS_PVT_WATER_H_
