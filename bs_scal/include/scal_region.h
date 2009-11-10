/** 
 * \file scal_region.h
 * \brief scal_region definition
 * \author Sergey Miryanov
 * \date 30.12.2008
 * */
#ifndef BS_SCAL_SCAL_REGION_H_
#define BS_SCAL_SCAL_REGION_H_

#include "scal_interpolate.h"
#include "scal_data_vector.h"

namespace blue_sky {

  template <typename strategy_t>
  class scale_array_holder;

  template <typename strategy_t>
  class scal_region
  {
  public:
    //typedef data_vector <strategy_t>          data_vector_t;
    typedef scale_array_holder <strategy_t>   scale_array_holder_t;
    typedef scal_region <strategy_t>          this_t;
    typedef scal_region_info <strategy_t>     scal_region_info_t;
    typedef typename strategy_t::item_t       item_t;

    scal_region (const scal_region_info_t &info_,
                 const data_vector <strategy_t> &Sp_,
                 const data_vector <strategy_t> &So_,
                 const data_vector <strategy_t> &Krp_,
                 const data_vector <strategy_t> &Krop_,
                 const data_vector <strategy_t> &Pcp_
                )
      : Sp (Sp_)
      , So (So_)
      , Krp (Krp_)
      , Krop (Krop_)
      , Pcp (Pcp_)
      , info (info_)
    {
    }

    /**
     * \brief do work
     *
     * \param[in] cell_index
     * \param[in] sat
     * \param[in] scale_arrays
     * \param[out] kr
     * \param[out] d_kr
     * \param[out] kro
     * \param[out] d_kro
     * */
    void process_2phase (int cell_index, item_t sat, item_t oil_sat, const scale_array_holder_t &scale_arrays,
                         item_t so_sub, item_t t_so_sub,
                         item_t &kr, item_t &d_kr, item_t &kro, item_t &d_kro) const
    {
      item_t s_max      = get_phase_sat_max ();
      item_t s_min      = get_phase_sat_min ();
      item_t so_cr      = get_sorp ();
      item_t s_cr       = get_spr ();

      item_t socr       = scale_arrays.get_socr (so_cr) [cell_index];
      item_t scr        = scale_arrays.get_scr (s_cr)   [cell_index];
      item_t su         = scale_arrays.get_su (s_max)   [cell_index];
      item_t sl         = scale_arrays.get_sl (s_min)   [cell_index];

      item_t t_so_max   = 1.0 - sl - t_so_sub;
      item_t so_max     = 1.0 - s_min - so_sub;

      item_t s          = scale_table (scr,  s_cr,  su,       s_max,  sat);
      item_t so         = scale_table (socr, so_cr, t_so_max, so_max, oil_sat);

      item_t d_kr_mult  = (s_max    - s_cr)  / (su       - scr);
      item_t d_kro_mult = (t_so_max - so_cr) / (t_so_max - socr);

      interpolate (s, Sp, Krp, s_cr, info.Krp_min_greater_zero, kr, d_kr, std::less <item_t> ());
      interpolate (so, So, Krop, so_cr, info.Krop_min_greater_zero, kro, d_kro, std::less <item_t> ());

      d_kr  = d_kr * d_kr_mult;
      d_kro = /*-*/d_kro * d_kro_mult;
    }

    void process_capillary (int cell_index, item_t sat, const scale_array_holder_t &scale_arrays,
                            item_t &cap, item_t *d_cap) const
    {
      item_t s_max  = get_phase_sat_max ();
      item_t s_min  = get_phase_sat_min ();

      item_t su     = scale_arrays.get_su (s_max) [cell_index];
      item_t sl     = scale_arrays.get_sl (s_min) [cell_index];

      item_t s      = scale_table (sl, s_min, su, s_max, sat);
      item_t mult   = (s_max - s_min) / (su - sl);

      interpolate (s, Sp, Pcp, cap, d_cap, std::less <item_t> ());

      if (d_cap)
        *d_cap = *d_cap * mult;
    }

    void process_init (int cell_index, const item_t cap, const scale_array_holder_t &scale_arrays, item_t &sat) const
    {
      item_t s;

      interpolate (Pcp, Sp, cap, s, std::less <item_t> ());

      item_t s_max  = get_phase_sat_max ();
      item_t s_min  = get_phase_sat_min ();

      item_t su     = scale_arrays.get_su (s_max) [cell_index];
      item_t sl     = scale_arrays.get_sl (s_min) [cell_index];

      sat = scale_not_table (sl, s_min, su, s_max, s);
    }

    item_t  get_spr () const
    {
      return info.spr;
    }
    item_t  get_sorp () const
    {
      return info.sorp;
    }
    item_t  get_kpr () const
    {
      return info.kpr;
    }
    item_t  get_krorp () const
    {
      return info.krorp;
    }
    item_t  get_phase_sat_max () const
    {
      return Sp.back ();
    }
    item_t  get_phase_sat_min () const
    {
      return Sp.front ();
    }
    item_t  get_pcp_max () const
    {
      return info.pcp_max;
    }

  public:

    data_vector <strategy_t> Sp;
    data_vector <strategy_t> So;
    data_vector <strategy_t> Krp;
    data_vector <strategy_t> Krop;
    data_vector <strategy_t> Pcp;

  private:
    const scal_region_info_t &info;
  };

}


#endif  // #ifndef BS_SCAL_SCAL_REGION_H_
