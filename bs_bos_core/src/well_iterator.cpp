/**
 * \file well_iterator.cpp
 * \brief impl of
 * \author Sergey Miryanov
 * \date 08.08.2008
 * */
#include "stdafx.h"
//
//#include "well_iterator.h"
//#include "facility_manager.h"
//#include "calc_well.h"
//#include "well_connection.h"
//#include "reservoir.h"
//#include "facility_manager.h"
//namespace blue_sky
//  {
//
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  well_iterator<list_t, iterator_t, value_type_t>::well_iterator (const list_t &list_, const iterator_t &it_)
//      : list_ (list_)
//      , it_ (it_)
//  {
//  }
//
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  bool
//  well_iterator<list_t, iterator_t, value_type_t>::operator == (const this_t &rhs) const
//    {
//      return it_ == rhs.it_;
//    }
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  bool
//  well_iterator<list_t, iterator_t, value_type_t>::operator != (const this_t &rhs) const
//    {
//      return it_ != rhs.it_;
//    }
//
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  const typename well_iterator<list_t, iterator_t, value_type_t>::this_t&
//  well_iterator<list_t, iterator_t, value_type_t>::operator = (const this_t &src)
//  {
//    list_->clear();
//    list_->insert (src.list_->begin(),src.list_->end());
//    it_ = src.it_;
//    return *this;
//  }
//
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  typename well_iterator<list_t, iterator_t, value_type_t>::value_type
//  well_iterator<list_t, iterator_t, value_type_t>::operator * () const
//    {
//      typename value_type_t::base_t w (it_->second, bs_dynamic_cast ());
//      BS_ASSERT (w);
//      if (!w)
//        throw bs_exception ("well_iterator::operator *()", "Oh, iterator value is not a well");
//
//      return w;
//    }
//
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  typename well_iterator<list_t, iterator_t, value_type_t>::this_t &
//  well_iterator<list_t, iterator_t, value_type_t>::operator ++ ()
//  {
//    iterator_t end = list_->end ();
//
//    ++it_;
//    while (it_ != end)
//      {
//        typename value_type_t::base_t w (it_->second, bs_dynamic_cast ());
//
//        if (w)
//          break;
//        else
//          ++it_;
//      }
//
//    return *this;
//  }
//
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  typename well_iterator<list_t, iterator_t, value_type_t>::this_t
//  well_iterator<list_t, iterator_t, value_type_t>::operator ++ (int)
//  {
//    this_t it = *this;
//
//    ++it;
//    return it;
//  }
//
//  //////////////////////////////////////////////////////////////////////////
//  typedef facility_manager <base_strategy_fi> fm_fi;
//  typedef facility_manager <base_strategy_di> fm_di;
//
//  template class well_iterator <fm_fi::locked_facility_map_t, fm_fi::facility_iterator_t, fm_fi::locked_well_t>;
//  template class well_iterator <fm_di::locked_facility_map_t, fm_di::facility_iterator_t, fm_di::locked_well_t>;
//
//} // namespace blue_sky
//
