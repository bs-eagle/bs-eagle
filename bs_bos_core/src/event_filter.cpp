/**
 *       \file  event_filter.cpp
 *      \brief  Impementation of event_filter
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  15.12.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "event_filter.h"
#include "bs_kernel.h"

namespace blue_sky
  {

  event_filter::event_filter(bs_type_ctor_param /*param = NULL */)
  {
  }
  event_filter::event_filter(const event_filter &x)
  : bs_refcounter (x), objbase (x)
  {

  }

  bool
  event_filter::accept_well (const std::string &well_name) const
    {
      if (reject_all_)
        return false;

      well_name_list_t::const_iterator it = std::find (well_name_list_.begin (), well_name_list_.end (), well_name);
      if (it != well_name_list_.end ())
        return false;

      return true;
    }

  void
  event_filter::add_filter_well (const std::string &well_name)
  {
    well_name_list_t::iterator it = std::find (well_name_list_.begin (), well_name_list_.end (), well_name);
    if (it == well_name_list_.end ())
      {
        well_name_list_.push_back (well_name);
      }
    else
      {
        BS_ASSERT (false && "well_name already in list") (well_name);
      }
  }

  void
  event_filter::set_reject_all (bool reject_all)
  {
    reject_all_ = reject_all;
  }

  BLUE_SKY_TYPE_STD_CREATE (event_filter);
  BLUE_SKY_TYPE_STD_COPY (event_filter);
  BLUE_SKY_TYPE_IMPL_SHORT (event_filter, objbase, "event_filter");

} // namespace blue_sky

