/**
 *       \file  default_connection.h
 *      \brief  Default implementation of well connection
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  20.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Should be moved to src/ directory
 * */
#ifndef BS_BOS_CORE_DEFAULT_CONNECTION_H_
#define BS_BOS_CORE_DEFAULT_CONNECTION_H_

#include "well_connection.h"
#include "shared_vector.h"

namespace blue_sky {
namespace wells {

  /**
   * \class default_connection
   * \brief Default implementation of well connection
   * */
  class BS_API_PLUGIN default_connection : public connection 
  {
  public:

    typedef connection base_t;
    typedef base_t::item_t             item_t;
    typedef base_t::rhs_item_t         rhs_item_t;

  public:

    //! blue-sky type declaration
    BLUE_SKY_TYPE_DECL (default_connection);

    /**
     * \brief  Clears data
     * */
    void 
    clear_data ();

    /**
     * \brief  Returns rw value
     * \return rw_value
     * \todo   Obsolete, should be removed
     * */
    shared_vector <item_t> 
    get_rw_value ();

    /**
     * \brief  Returns wr value
     * \return wr_value
     * \todo   Obsolete, should be removed
     * */
    shared_vector <item_t> 
    get_wr_value ();

    /**
     * \brief  Returns rr value
     * \return rr_value
     * \todo   Obsolete, should be removed
     * */
    shared_vector <item_t> 
    get_rr_value ();

    /**
     * \brief  Returns ps value
     * \return ps_value
     * \todo   Obsolete, should be removed
     * */
    shared_vector <item_t> 
    get_ps_value ();

    /**
     * \brief  Returns rate value
     * \return rate_value
     * \todo   Obsolete, should be removed
     * */
    shared_vector <rhs_item_t> 
    get_rate_value ();

  public:

    enum 
      {
        rr_value_count = FI_PHASE_TOT * FI_PHASE_TOT + FI_PHASE_TOT,      //!< rr_value now stores rr values and ps values
      };

    boost::array <item_t, FI_PHASE_TOT>                 mobility_value;   //!< Mobility value
    boost::array <rhs_item_t, FI_PHASE_TOT>             rate_value;       //!< Rate value
    boost::array <item_t, FI_PHASE_TOT>                 ps_value;         //!< \todo Obsolete, should be removed
    boost::array <item_t, rr_value_count>               rr_value;         //!< RR and PS values
    boost::array <item_t, FI_PHASE_TOT>                 rw_value;         //!< RW value
    boost::array <item_t, FI_PHASE_TOT>                 wr_value;         //!< WR value
  };


} // namespace wells
} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_DEFAULT_CONNECTION_H_


