/**
 * \file scale_array_holder.h
 * \brief
 * \author Sergey Miryanov
 * \date 02.12.2008
 * */
#ifndef BS_SCAL_SCALE_ARRAY_HOLDER_H_
#define BS_SCAL_SCALE_ARRAY_HOLDER_H_

#include "scale_arrays_placement_strategies.h"
#include "array_ext.h"

namespace blue_sky
  {

#ifdef _DEBUG
  struct value_accessor
    {
      value_accessor (const float *array_, int step_, int count_, double value_)
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

          return step_ ? array_[index * step_] : value_;
        }

private:

      const float   *array_;
      int           step_;
      int           count_;
      double        value_;

    };
#else
  struct value_accessor
    {
      value_accessor (const float *array__, int step__, int, double value__)
          : array_ (array__)
          , step_ (step__)
          , value_ (value__)
      {

      }

      double operator[] (int index) const
        {
          return step_ ? array_[index * step_] : value_;
        }

private:
      const float    *array_;
      int             step_;
      double          value_;
    };
#endif

  template <typename strategy_t>
  class BS_API_PLUGIN scale_array_holder : public objbase
    {
    public:
      typedef scale_array_holder<strategy_t>  this_t;

      typedef typename strategy_t::item_t     item_t;
      typedef typename strategy_t::item_array_t vector_t;
      typedef void (this_t::*inserter_t) (const array_ext <float> &);
      typedef typename strategy_t::template vec <float> vec_t;
      typedef typename vec_t::type             data_t;

    public:
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

    private:

      typedef scal::data_placement::separate_scale_vectors_t placement_t;
      //typedef scal::data_placement::struct_like_scale_arrays_t placement_t;

      data_t data;               ///< store all data
      scal::data_placement::scale_array_placement_info placement_info;
    public:

      BLUE_SKY_TYPE_DECL_T (scale_array_holder);
    };


} // namespace blue_sky



#endif  // #ifndef BS_SCAL_SCALE_ARRAY_HOLDER_H_

