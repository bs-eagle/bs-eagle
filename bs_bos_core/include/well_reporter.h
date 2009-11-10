/**
 * \file well_reporter.h
 * \brief
 * \date 01.10.2009
 * \author Sergey Miryanov
 * */
#ifndef BS_BOS_CORE_WELL_REPORTER_H_
#define BS_BOS_CORE_WELL_REPORTER_H_

#include BS_FORCE_PLUGIN_IMPORT ()
#include "data_class.h"
#include BS_STOP_PLUGIN_IMPORT ()
#include "reservoir.h"

namespace blue_sky {

  template <typename strategy_t>
  struct BS_API_PLUGIN well_data_printer
  {
    static void
    print_prod (const smart_ptr <idata, true> &data, const smart_ptr <reservoir <strategy_t>, true> &rs);
    static void
    print_inj  (const smart_ptr <idata, true> &data, const smart_ptr <reservoir <strategy_t>, true> &rs);
    static void
    print_total_prod (const smart_ptr <idata, true> &data, const smart_ptr <reservoir <strategy_t>, true> &rs);
  };
  

} // namespace blue_sky


#endif // #ifndef BS_BOS_CORE_WELL_REPORTER_H_

