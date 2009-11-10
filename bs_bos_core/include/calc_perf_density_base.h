/**
 * \file calc_perf_density_base.h
 * \brief base class for perforation's density calculation
 * \author Sergey Miryanov
 * \date 18.11.2008
 * */
#ifndef BS_CALC_PERF_DENSITY_BASE_H_
#define BS_CALC_PERF_DENSITY_BASE_H_


namespace blue_sky
  {

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  class well;

  template <typename strategy_t>
  class calc_perf_density_base : public objbase
    {
    public:

      typedef typename strategy_t::index_t      index_t;
      typedef typename strategy_t::item_t       item_t;

      typedef calc_model <strategy_t>           calc_model_t;
      typedef well <strategy_t>                 well_t;

      typedef smart_ptr <calc_model_t, true>    sp_calc_model_t;
      typedef smart_ptr <well_t, true>          sp_well_t;

    public:
      virtual ~calc_perf_density_base () {}

      virtual void
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const = 0;
    };


} // namespace blue_sky


#endif // #ifndef BS_CALC_PERF_DENSITY_BASE_H_

