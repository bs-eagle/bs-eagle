/**
 *       \file  well_reporter.h
 *      \brief  Functions to print well information tables
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  01.10.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_BOS_CORE_WELL_REPORTER_H_
#define BS_BOS_CORE_WELL_REPORTER_H_

#include "bs_common.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "data_class.h"
#include BS_STOP_PLUGIN_IMPORT ()
#include "reservoir.h"

namespace blue_sky {

  /**
   * \class well_data_printer
   * \brief Prints well information tables
   * */
  struct BS_API_PLUGIN well_data_printer
  {
    /**
     * \brief  Prints production information for all wells
     * \param  data Data holder
     * \param  rs Instance of reservoir_simulator
     * */
    static void
    print_prod (const smart_ptr <idata, true> &data, const smart_ptr <reservoir, true> &rs);

    /**
     * \brief  Prints injection information for all wells
     * \param  data Data holder
     * \param  rs Instance of reservoir_simulator
     * */
    static void
    print_inj  (const smart_ptr <idata, true> &data, const smart_ptr <reservoir, true> &rs);

    /**
     * \brief  Prints total production and injection information
     *         for all wells
     * \param  data Data holder
     * \param  rs Instance of reservoir_simulator
     * \todo   Rename
     * */
    static void
    print_total_prod (const smart_ptr <idata, true> &data, const smart_ptr <reservoir, true> &rs);
  };
  

} // namespace blue_sky


#endif // #ifndef BS_BOS_CORE_WELL_REPORTER_H_

