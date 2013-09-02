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
  class BS_API_PLUGIN pvt_gas : public pvt_base 
    {
    public:

      typedef pvt_base                 base_t;
      typedef base_t::vector_t         vector_t;

      enum {
         PVT_GAS_INPUT_PRESSURE = 0,
         PVT_GAS_INPUT_FVF,
         PVT_GAS_INPUT_VISC,
         PVT_GAS_INPUT_TOTAL
       };  

      enum {
         PVT_GAS_PRESSURE = 0,
         PVT_GAS_INV_FVF,
         PVT_GAS_INV_VISC,
         PVT_GAS_INV_VISC_FVF,
         PVT_GAS_TOTAL
      };

      /**
       * \brief store values into data
       *
       * \param seq_vector
       */
      virtual void insert_vector (const v_double &vec);

      /**
       * \brief generate interpolated data
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
      virtual void calc (const t_double p, t_float  *inv_fvf, t_float *d_inv_fvf,
                         t_float *inv_visc, t_float *d_inv_visc,
                         t_float *inv_visc_fvf, t_float *d_inv_visc_fvf) const;

      virtual void
      print () const;

    private:

      void check_gas ();
      void check_gas_internal ();

    private:

      //vector_t main_pressure_;
      //vector_t main_fvf_;
      //vector_t main_visc_;

      //vector_t pressure_;
      //vector_t inv_fvf_;
      //vector_t inv_visc_;
      //vector_t inv_visc_fvf_;

    public:

      BLUE_SKY_TYPE_DECL_T (pvt_gas);
    };

} // namespace blue_sky


#endif  // #ifndef BS_PVT_GAS_H_
