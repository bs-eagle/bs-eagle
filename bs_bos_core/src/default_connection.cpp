/**
 *       \file  default_connection.cpp
 *      \brief  Implementation of default connection
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  20.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "default_connection.h"
#include "vector_assign.h"

namespace blue_sky {
namespace wells {

  /**
   * \brief  'default' ctor for default_connection
   * \param  param Additional params for ctor
   * */
  default_connection::default_connection (bs_type_ctor_param param /* = NULL */)
  : base_t (param)
  {
    clear_data ();
  }
  /**
   * \brief  copy-ctor for default_connection
   * \param  c Instance of default_connection to be copied
   * */
  default_connection::default_connection (const default_connection &c)
        : bs_refcounter (), connection() 
  {
    *this = c;
  }

  shared_vector <default_connection::item_t>
  default_connection::get_rw_value ()
  {
    return shared_array <item_t> (&rw_value[0], rw_value.size ());
  }
  shared_vector <default_connection::item_t>
  default_connection::get_wr_value ()
  {
    return shared_array <item_t> (&wr_value[0], wr_value.size ());
  }
  shared_vector <default_connection::item_t>
  default_connection::get_rr_value ()
  {
    return shared_array <item_t> (&rr_value[0], rr_value.size ());
  }
  shared_vector <default_connection::item_t>
  default_connection::get_ps_value ()
  {
    return shared_array <item_t> (&ps_value[0], ps_value.size ());
  }
  shared_vector <default_connection::rhs_item_t>
  default_connection::get_rate_value ()
  {
    return shared_array <rhs_item_t> (&rate_value[0], rate_value.size ());
  }

  void
  default_connection::clear_data ()
  {
    assign (mobility_value, 0);
    assign (rate_value, 0);
    assign (ps_value, 0);
    assign (rr_value, 0);
    assign (rw_value, 0);
    assign (wr_value, 0);

    base_t::clear_data ();
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (default_connection);
  BLUE_SKY_TYPE_STD_COPY (default_connection);
  BLUE_SKY_TYPE_IMPL (default_connection, connection, "default_connection_seq", "default_connection_seq", "default_connection_seq");
  //////////////////////////////////////////////////////////////////////////

} //namespace wells
} //namespace blue_sky

