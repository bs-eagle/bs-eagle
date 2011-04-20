/**
 * \file scale_array_holder.h
 * \brief
 * \author Sergey Miryanov
 * \date 02.12.2008
 * */
#ifndef BS_SCAL_SCALE_ARRAY_HOLDER_H_
#define BS_SCAL_SCALE_ARRAY_HOLDER_H_

#include "scale_arrays_placement_strategies.h"
#include "shared_vector.h"

namespace blue_sky
  {

  struct value_accessor
    {
      value_accessor (const t_float *array_, int step_, int count_, double value_)
          : array_ (array_)
          , step_ (step_)
          , count_ (count_)
          , value_ (value_)
      {

      }

      double operator [] (int index) const
        {
#ifdef _DEBUG
          if (step_)
            {
              BS_ASSERT (index >= 0 && index < count_);
            }
#endif

          return count_ ? array_[index * step_] : value_;
        }

private:

      const t_float   *array_;
      int           step_;
      int           count_;
      double        value_;

    };

  class BS_API_PLUGIN scale_array_holder : public objbase
    {
    public:
      typedef scale_array_holder                this_t;

      typedef t_double                          item_t;
      //typedef v_double                        vector_t;
      typedef void (this_t::*inserter_t) (const shared_vector <float> &);
      typedef shared_vector <t_float>           data_t;

      typedef vartype_table_iface <t_float>              vartype_table_iface_t;
      typedef vartype_table_iface_t::vector_t            vector_t;
      typedef smart_ptr<vartype_table_iface_t,true>      sp_vartype_table; 

      enum array_name
        {
          socr,
          scr,
          su,
          sl,
          pcp,
          krp,
          krop,
          krpr,
          krorp,
          keyword_total
        };

    public:

#if 0    
      ///< critical oil saturation (NX * NY * NZ)
      value_accessor get_socr (item_t value_) const
      {
        return value_accessor (&data[placement_info.socr_offset], placement_info.socr_step, placement_info.size, value_);
      }
      ///< critical (water/gas) saturation (NX * NY * NZ)
      value_accessor get_scr  (item_t value_) const
      {
        return value_accessor (&data[placement_info.scr_offset], placement_info.scr_step, placement_info.size, value_);
      }
      ///< max saturation value in SCAL table (NX * NY * NZ)
      value_accessor get_su   (item_t value_) const
      {
        return value_accessor (&data[placement_info.su_offset], placement_info.su_step, placement_info.size, value_);
      }
      ///< min saturation value in SCAL table (NX * NY * NZ)
      value_accessor get_sl   (item_t value_) const
      {
        return value_accessor (&data[placement_info.sl_offset], placement_info.sl_step, placement_info.size, value_);
      }
      ///< max capillary pressure value in SCAL table for water oil system (NX * NY * NZ), for water only
      value_accessor get_pcp  (item_t value_) const
      {
        return value_accessor (&data[placement_info.pcp_offset], placement_info.pcp_step, placement_info.size, value_);
      }
#else
      ///< critical oil saturation (NX * NY * NZ)
      value_accessor get_socr (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (socr), 1, data_pool->get_n_rows (socr), value_);
      }
      ///< critical (water/gas) saturation (NX * NY * NZ)
      value_accessor get_scr  (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (scr), 1, data_pool->get_n_rows (scr), value_);
      }
      ///< max saturation value in SCAL table (NX * NY * NZ)
      value_accessor get_su   (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (su), 1, data_pool->get_n_rows (su), value_);
      }
      ///< min saturation value in SCAL table (NX * NY * NZ)
      value_accessor get_sl   (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (sl), 1, data_pool->get_n_rows (sl), value_);
      }
      ///< max capillary pressure value in SCAL table for water oil system (NX * NY * NZ), for water only
      value_accessor get_pcp  (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (pcp), 1, data_pool->get_n_rows (pcp), value_);
      }
      ///< max scaled end point phase relative permeability in SCAL table for phase-oil system (NX * NY * NZ)
      value_accessor get_krp  (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (krp), 1, data_pool->get_n_rows (krp), value_);
      }
      ///< max scaled end point oil relative permeability in SCAL table for phase-oil system (NX * NY * NZ)
      value_accessor get_krop  (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (krop), 1, data_pool->get_n_rows (krop), value_);
      }

      ///< scaled end point phase relative permeability at residual oil saturation in SCAL table for phase-oil system (NX * NY * NZ)
      value_accessor get_krpr  (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (krpr), 1, data_pool->get_n_rows (krpr), value_);
      }
      
      ///< scaled end point phase relative permeability at critical phase saturation in SCAL table for phase-oil system (NX * NY * NZ)
      value_accessor get_krorp  (item_t value_) const
      {
        return value_accessor (data_pool->get_col_ptr (krorp), 1, data_pool->get_n_rows (krorp), value_);
      }

#endif 
      template <typename array_t>
      void insert_socr (const array_t &socr_)
      {
        if (socr_.size ())
          placement_t::place_data (placement_t::socr, data, socr_, placement_info);
      }
      template <typename array_t>
      void insert_scr (const array_t &scr_)
      {
        if (scr_.size ())
          placement_t::place_data (placement_t::scr, data, scr_, placement_info);
      }
      template <typename array_t>
      void insert_su (const array_t &su_)
      {
        if (su_.size ())
          placement_t::place_data (placement_t::su, data, su_, placement_info);
      }
      template <typename array_t>
      void insert_sl (const array_t &sl_)
      {
        if (sl_.size ())
          placement_t::place_data (placement_t::sl, data, sl_, placement_info);
      }
      template <typename array_t>
      void insert_pcp (const array_t &pcp_)
      {
        if (pcp_.size ())
          placement_t::place_data (placement_t::pcp, data, pcp_, placement_info);
      }
      void remove_pcp ()
      {
        placement_t::remove_data (placement_t::pcp, data, placement_info);
      }
      
      template <typename array_t>
      void set_socr (const array_t &socr_)
        {
          if (socr_->size ())
            data_pool->add_col_vector (socr, "SOCR", socr_);
        }

      template <typename array_t>
      void set_scr (const array_t &scr_)
        {
          if (scr_->size ())
            data_pool->add_col_vector (scr, "SCR", scr_);
        }

      template <typename array_t>
      void set_su (const array_t &su_)
        {
          if (su_->size ())
            data_pool->add_col_vector (su, "SU", su_);
        }

      template <typename array_t>
      void set_sl (const array_t &sl_)
        {
          if (sl_->size ())
            data_pool->add_col_vector (sl, "SL", sl_);
        }

      template <typename array_t>
      void set_pcp (const array_t &pcp_)
        {
          if (pcp_->size ())
            data_pool->add_col_vector (pcp, "PCP", pcp_);
        }

      template <typename array_t>
      void set_krp (const array_t &krp_)
        {
          if (krp_->size ())
            data_pool->add_col_vector (krp, "KRP", krp_);
        }

      template <typename array_t>
      void set_krop (const array_t &krop_)
        {
          if (krop_->size ())
            data_pool->add_col_vector (krop, "KROP", krop_);
        }

      template <typename array_t>
      void set_krpr (const array_t &krpr_)
        {
          if (krpr_->size ())
            data_pool->add_col_vector (krpr, "KRPR", krpr_);
        }

      template <typename array_t>
      void set_krorp (const array_t &krorp_)
        {
          if (krorp_->size ())
            data_pool->add_col_vector (krorp, "KRORP", krorp_);
        }
 
      int 
      is_prop_valid (const t_long col) const
        {
          return (data_pool->get_n_rows (col) == 0) ? 0 : 1;
        }
      
    private:

      typedef scal::data_placement::separate_scale_vectors_t placement_t;
      //typedef scal::data_placement::struct_like_scale_arrays_t placement_t;

      data_t data;               ///< store all data
      scal::data_placement::scale_array_placement_info placement_info;
      sp_vartype_table data_pool;
    public:

      BLUE_SKY_TYPE_DECL_T (scale_array_holder);
    };


} // namespace blue_sky



#endif  // #ifndef BS_SCAL_SCALE_ARRAY_HOLDER_H_

