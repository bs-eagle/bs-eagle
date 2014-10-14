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
#include "scal_region_info.h"

namespace blue_sky {

  class scale_array_holder;

  class scal_region
  {
  public:
    //typedef data_vector <strategy_t>          data_vector_t;
    typedef t_float                           item_t;
    typedef scale_array_holder                scale_array_holder_t;
    typedef scal_region                       this_t;
    typedef scal_region_info <item_t>         scal_region_info_t;

    scal_region (const scal_region_info_t &info_,
                 const data_vector <item_t> &Sp_,
                 const data_vector <item_t> &So_,
                 const data_vector <item_t> &Krp_,
                 const data_vector <item_t> &Krop_,
                 const data_vector <item_t> &Pcp_
                )
      : Sp (Sp_)
      , So (So_)
      , Krp (Krp_)
      , Krop (Krop_)
      , Pcp (Pcp_)
      , info (info_)
    {
    }


    void process_vertical_scaling (t_long /*cell_index*/, item_t sat, item_t sat_scale,
                                   item_t su, item_t t_sr_phase,
                                   item_t t_krp, item_t krp_max,
                                   item_t t_krpr, item_t krp_sr,
                                   int is_krpr,
                                   item_t &kr, item_t &d_kr)
    {
      if (is_krpr)
        {
          item_t kr_mult = 1.0;
            
          if (sat < t_sr_phase)
            {
              kr_mult = t_krpr / krp_sr;
              kr = kr * kr_mult;
            }   
          else
            {
              if (fabs (krp_max - krp_sr) > EPS_DIV)
                {
                  kr_mult = (t_krp - t_krpr) / (krp_max - krp_sr);
                  kr = t_krpr + (kr - krp_sr) * kr_mult;
                }
              else //  fabs (krp_max - krp_sr) < EPS_DIV 
                {
                  // TODO: make linear function between KRP and KRPR
                  if (fabs (su - t_sr_phase) > EPS_DIV)
                    {
                      kr_mult = (t_krp - t_krpr) / (su - t_sr_phase);
                      kr = t_krpr + (sat_scale - t_sr_phase) * kr_mult;
                    }
                  else // 
                    {
                       // TODO :
                       BS_ASSERT (false && "KRPR scaling case (SPU = SR): NOT IMPL YET");
                    }    
                }
            } 
          d_kr = d_kr * kr_mult;   
        }
      else 
        {
          item_t kr_mult = t_krp / krp_max;
          kr = kr * kr_mult;
          d_kr = d_kr * kr_mult;
        }   
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
    void process_2phase (t_long cell_index, item_t sat, item_t oil_sat, const scale_array_holder_t &scale_arrays,
                         item_t so_sub, item_t t_so_sub,
                         item_t &kr, item_t &d_kr, item_t &kro, item_t &d_kro,
                         int is_crs) const
    {
      item_t s_max      = get_phase_sat_max ();
      item_t s_min      = get_phase_sat_min ();
      item_t so_cr      = get_sorp ();
      item_t s_cr       = get_spr ();

      item_t socr       = scale_arrays.get (blue_sky::socr, so_cr) [cell_index];
      item_t scr        = scale_arrays.get (blue_sky::scr, s_cr)   [cell_index];
      item_t su         = scale_arrays.get (blue_sky::su, s_max)   [cell_index];
      item_t sl         = scale_arrays.get (blue_sky::sl, s_min)   [cell_index];

      item_t t_so_max   = 1.0 - sl - t_so_sub;
      item_t so_max     = 1.0 - s_min - so_sub;

      item_t sr_phase   = 1.0 - so_cr - so_sub;  
      item_t t_sr_phase = 1.0 - socr - t_so_sub;
      
      item_t sr_oil     = 1.0 - s_cr - so_sub;
      item_t t_sr_oil   = 1.0 - scr - t_so_sub;
      
      item_t s          = 0.0;
      item_t so         = 0.0;
      item_t d_kr_mult  = 0.0;
      item_t d_kro_mult = 0.0;
      
      if (is_crs)
        {
          // phase saturation scale
          if (sat < t_sr_phase)
            {
              s         = scale_table (scr, s_cr, t_sr_phase, sr_phase,  sat);
              d_kr_mult = (sr_phase - s_cr) / (t_sr_phase - scr);
            }
          else
            {
              s         = scale_table (t_sr_phase, sr_phase, su, s_max,  sat);       
              d_kr_mult = (s_max - sr_phase) / (su - t_sr_phase);     
            }  
          // oil saturation scale
          if (oil_sat < t_sr_oil)
            {
              so         = scale_table (socr, so_cr, t_sr_oil, sr_oil, oil_sat);
              d_kro_mult = (sr_oil - so_cr) / (t_sr_oil - socr);
            }  
          else
            {
              so         = scale_table (t_sr_oil, sr_oil, t_so_max, so_max,  oil_sat); 
              d_kro_mult = (so_max - sr_oil) / (su - t_sr_phase);
            }  
        }      
      else
        {
          s          = scale_table (scr,  s_cr,  su,       s_max,  sat);
          so         = scale_table (socr, so_cr, t_so_max, so_max, oil_sat);
          d_kr_mult  = (s_max    - s_cr)  / (su       - scr);
          d_kro_mult = (t_so_max - so_cr) / (t_so_max - socr);
        }  
      
      interpolate (s, Sp, Krp, s_cr, info.Krp_min_greater_zero, kr, d_kr, std::less <item_t> ());
      interpolate (so, So, Krop, so_cr, info.Krop_min_greater_zero, kro, d_kro, std::less <item_t> ());

      d_kr  = d_kr * d_kr_mult;
      d_kro = /*-*/d_kro * d_kro_mult;

      ////////////////////////////////////////////////////////////
      // Phase relative permeability (vertical) scaling
      ////////////////////////////////////////////////////////////
      {
        item_t krp_max = get_krp_max ();
        item_t t_krp = scale_arrays.get (blue_sky::krp, krp_max) [cell_index];
      
        if (scale_arrays.is_prop_valid (blue_sky::krpr))
          {
            item_t kr_mult = 1.0;
            item_t krp_sr;  
            item_t t_krpr;
            if (is_crs)
              {
                (void)t_sr_phase;
                //t_sr_phase = t_sr_phase;
              }
            else
              {
                t_sr_phase = scale_table (scr, s_cr, su, s_max, sr_phase);
              }  
              
            krp_sr = calc_krp (t_sr_phase);
            t_krpr = scale_arrays.get (blue_sky::krpr, krp_sr) [cell_index]; 
            if (sat < t_sr_phase)
              {
                kr_mult = t_krpr / krp_sr;
                kr = kr * kr_mult;
              }   
            else
              {
                if (fabs (krp_max - krp_sr) > EPS_DIV)
                  {
                    kr_mult = (t_krp - t_krpr) / (krp_max - krp_sr);
                    kr = t_krpr + (kr - krp_sr) * kr_mult;
                  }
                else //  fabs (krp_max - krp_sr) < EPS_DIV 
                  {
                    // TODO: make linear function between KRP and KRPR
                    if (fabs (su - t_sr_phase) > EPS_DIV)
                      {
                        kr_mult = (t_krp - t_krpr) / (su - t_sr_phase);
                        kr = t_krpr + (s - t_sr_phase) * kr_mult;
                      }
                    else // 
                      {
                         // TODO :
                         BS_ASSERT (false && "KRPR scaling case (SPU = SR): NOT IMPL YET");
                      }    
                  }
              } 
            d_kr = d_kr * kr_mult;   
          }
        else 
          {
            item_t kr_mult = t_krp / krp_max;
            kr = kr * kr_mult;
            d_kr = d_kr * kr_mult;
          }   
      }    
      ////////////////////////////////////////////////////////////
      // Oil relative permeability (vertical) scaling
      ////////////////////////////////////////////////////////////
      {
        item_t krop_max = get_krop_max ();
        item_t t_krop = scale_arrays.get (blue_sky::krop, krop_max) [cell_index];
        
        if (scale_arrays.is_prop_valid (blue_sky::krorp))
          {
            item_t krop_mult = 1.0;
            item_t krop_sr;  
            item_t t_krorp;
            if (is_crs)
              {
                (void)t_sr_oil;
                //t_sr_oil = t_sr_oil;
              }
            else
              {
                t_sr_oil = scale_table (socr, so_cr, t_so_max, so_max, sr_oil);
              }  
              
            krop_sr = calc_krorp (t_sr_oil);
            t_krorp = scale_arrays.get (blue_sky::krorp, krop_sr) [cell_index]; 
            if (oil_sat < t_sr_oil)
              {
                krop_mult = t_krorp / krop_sr;
                kro = kro * krop_mult;
              }   
            else
              {
                if (fabs (krop_max - krop_sr) > EPS_DIV)
                  {
                    krop_mult = (t_krop - t_krorp) / (krop_max - krop_sr);
                    kro = t_krorp + (kro - krop_sr) * krop_mult;
                  }
                else //  fabs (krop_max - krop_sr) < EPS_DIV 
                  {
                    // TODO: make linear function between KRP and KRPR
                    if (fabs (t_so_max - t_sr_oil) > EPS_DIV)
                      {
                        krop_mult = (t_krop - t_krorp) / (t_so_max - t_sr_oil);
                        kro = t_krorp + (so - t_sr_oil) * krop_mult;
                      }
                    else // 
                      {
                         // TODO :
                         BS_ASSERT (false && "KRORP scaling case (SO_MAX = SR): NOT IMPL YET");
                      }    
                  }
              } 
            d_kro = d_kro * krop_mult;   
          }
        else 
          {
            item_t krop_mult = t_krop / krop_max;
            kro = kro * krop_mult;
            d_kro = d_kro * krop_mult;
          }   
      }
    }

    void process_capillary (int cell_index, item_t sat, const scale_array_holder_t &scale_arrays,
                            item_t &cap, item_t *d_cap) const
    {
      item_t s_max  = get_phase_sat_max ();
      item_t s_min  = get_phase_sat_min ();

      const value_accessor &su_ = scale_arrays.get (blue_sky::su, s_max);
      const value_accessor &sl_ = scale_arrays.get (blue_sky::sl, s_min);

      item_t su     = su_ [cell_index];
      item_t sl     = sl_ [cell_index];

      item_t s      = scale_table (sl, s_min, su, s_max, sat);
      item_t mult   = (s_max - s_min) / (su - sl);

      interpolate (s, Sp, Pcp, cap, d_cap, std::less <item_t> ());

      if (d_cap)
        *d_cap = *d_cap * mult;
    }

    void process_capillary_2 (item_t sat, item_t &cap) const
    {
      interpolate (sat, Sp, Pcp, cap, (item_t *) NULL, std::less <item_t> ());
    }

    void process_init (int cell_index, const item_t cap, const scale_array_holder_t &scale_arrays, item_t &sat) const
    {
      item_t s;

      interpolate (Pcp, Sp, cap, s, std::less <item_t> ());

      item_t s_max  = get_phase_sat_max ();
      item_t s_min  = get_phase_sat_min ();

      item_t su     = scale_arrays.get (blue_sky::su, s_max) [cell_index];
      item_t sl     = scale_arrays.get (blue_sky::sl, s_min) [cell_index];

      sat = scale_not_table (sl, s_min, su, s_max, s, true);
    }

    void process_init_2 (const item_t cap, item_t &sat) const
    {
      interpolate (Pcp, Sp, cap, sat, std::less <item_t> ());
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
    item_t get_krp_max () const 
    {
      return Krp.back ();
    }
    item_t get_krop_max () const 
    {
      return Krop.back ();
    } 
    item_t  get_pcp_max () const
    {
      return info.pcp_max;
    }

    item_t calc_krp (item_t sr_phase) const
    {
      item_t kpr_;
      interpolate (Sp, Krp, sr_phase, kpr_, std::less <item_t> ());
      return kpr_; 
    }

    item_t calc_krorp (item_t sr_oil) const
    {
      item_t krorp_;
      interpolate (So, Krop, sr_oil, krorp_, std::less <item_t> ());     
      return krorp_; 
    }
    
  public:

    data_vector <item_t> Sp;
    data_vector <item_t> So;
    data_vector <item_t> Krp;
    data_vector <item_t> Krop;
    data_vector <item_t> Pcp;

  private:
    const scal_region_info_t &info;
  };

}


#endif  // #ifndef BS_SCAL_SCAL_REGION_H_
