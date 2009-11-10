/**
 * \file for_each_well.h
 * \brief well traversers
 * \author Sergey Miryanov
 * \date 13.02.2009
 * */
#ifndef BS_FOR_EACH_WELL_H_
#define BS_FOR_EACH_WELL_H_

namespace blue_sky {

  template <typename strategy_t, typename method_t>
  void
  for_each_well (facility_manager <strategy_t> &facility_list, method_t method)
  {
    typename facility_manager <strategy_t>::well_const_iterator_t wb = facility_list.wells_begin ();
    typename facility_manager <strategy_t>::well_const_iterator_t we = facility_list.wells_end ();
    typedef well <strategy_t> well_t;
    typedef smart_ptr <well_t, true> sp_well_t;

    for (; wb != we; ++wb)
      {
        sp_well_t well (wb->second, bs_dynamic_cast ());
        BS_ASSERT (well);

        if (well->get_connections_count ())
          {
            method (*well);
          }
      }
  }

  template <typename strategy_t, typename method_t>
  void
  for_each_facility (facility_manager <strategy_t> &facility_list, method_t method)
  {
    typename facility_manager <strategy_t>::well_const_iterator_t wb = facility_list.wells_begin ();
    typename facility_manager <strategy_t>::well_const_iterator_t we = facility_list.wells_end ();

    for (; wb != we; ++wb)
      {
        method (*wb->second);

        //if (well->get_connections_count ())
        //  {
        //    method (*well);
        //  }
      }
  }

  template <typename strategy_t, typename method_t>
  bool
  for_each_well_while_cond (facility_manager <strategy_t> &facility_list, method_t method)
  {
    typename facility_manager <strategy_t>::well_const_iterator_t wb = facility_list.wells_begin ();
    typename facility_manager <strategy_t>::well_const_iterator_t we = facility_list.wells_end ();
    typedef well <strategy_t> well_t;
    typedef smart_ptr <well_t, true> sp_well_t;

    for (; wb != we; ++wb)
      {
        sp_well_t well (wb->second, bs_dynamic_cast ());
        BS_ASSERT (well);
        if (well->get_connections_count () && method (*well))
          {
            return true;
          }
      }

    return false;
  }
  template <typename strategy_t, typename method_t>
  size_t
  for_each_well_acc (facility_manager <strategy_t> &facility_list, method_t method)
  {
    typename facility_manager <strategy_t>::well_const_iterator_t wb = facility_list.wells_begin ();
    typename facility_manager <strategy_t>::well_const_iterator_t we = facility_list.wells_end ();
    typedef well <strategy_t> well_t;
    typedef smart_ptr <well_t, true> sp_well_t;

    size_t acc = 0;
    for (; wb != we; ++wb)
      {
        sp_well_t well (wb->second, bs_dynamic_cast ());
        BS_ASSERT (well);

        if (well->get_connections_count ())
          {
            acc += method (*well);
          }
      }

    return acc;
  }

} // namespace blue_sky

#endif  // #ifndef BS_FOR_EACH_WELL_H_

