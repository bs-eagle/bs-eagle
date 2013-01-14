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
  class BS_API_PLUGIN pvt_dead_oil : public pvt_base 
    {
    public:

      typedef pvt_base                          base_t;
      typedef base_t::vector_t                  vector_t;
      
      enum {
         PVT_OIL_INPUT_GPR = 0,
         PVT_OIL_INPUT_PRESSURE,
         PVT_OIL_INPUT_FVF,
         PVT_OIL_INPUT_VISC,
         PVT_OIL_INPUT_TOTAL
       };  
           
      enum {
         PVT_OIL_PRESSURE = 0,
         PVT_OIL_INV_FVF,
         PVT_OIL_INV_VISC,
         PVT_OIL_INV_VISC_FVF,
         PVT_OIL_GOR,
         PVT_OIL_TOTAL
      };
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
            t_float *inv_fvf, t_float *d_inv_fvf, t_float *gor_d_inv_fvf,
            t_float *inv_visc, t_float *d_inv_visc, t_float *gor_d_inv_visc,
            t_float *inv_visc_fvf, t_float *d_inv_visc_fvf, t_float *gor_d_inv_visc_fvf,
            t_float *gas_oil_ratio, t_float *d_gas_oil_ratio,
            const t_double drsdt = -1.0, const t_double dt = 0,
            const t_double old_gas_oil_ratio = 0) const;

      virtual
      t_double get_gor_for_pressure (t_double pressure_data) const;

      virtual
      t_double interpolate_and_fix (t_double cell_pbub) const;

      const vector_t &
      get_pressure () const 
      {
        return pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
      }
      const vector_t &
      get_gor () const 
      {
        return pvt_props_table->get_col_vector (PVT_OIL_GOR);
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
      int build_internal (t_double atm_p, t_double min_p, t_double max_p, t_long n_intervals, bool is_pvto);

      /**
       * \brief calc saturated oil
       */
      bool
      calc_saturated_oil (const bool is_g, const int main_var, const t_double p, const t_double gor,
                          t_float *inv_fvf, t_float *d_inv_fvf, t_float *gor_d_inv_fvf,
                          t_float *inv_visc, t_float *d_inv_visc, t_float *gor_d_inv_visc,
                          t_float *inv_visc_fvf, t_float *d_inv_visc_fvf, t_float *gor_d_inv_visc_fvf,
                          t_float *gas_oil_ratio, t_float *d_gas_oil_ratio,
                          const t_double drsdt = -1.0, const t_double dt = 0,
                          const t_double old_gas_oil_ratio = 0) const;

    protected:


//      vector_t  main_gpr_;
//      vector_t  main_pressure_;
//      vector_t  main_fvf_;
//      vector_t  main_visc_;

//      vector_t  pressure_;
//      vector_t  inv_fvf_;
//      vector_t  inv_visc_;
//      vector_t  inv_visc_fvf_;
//      vector_t  gor_;

    public:

      BLUE_SKY_TYPE_DECL_T (pvt_dead_oil);

    };

} // namespace blue_sky

#endif  // #ifndef BS_PVT_DEAD_OIL_H_
