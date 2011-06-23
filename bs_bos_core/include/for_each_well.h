/**
 *       \file  for_each_well.h
 *      \brief  Calls methods for each facility (or well) in list
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  13.02.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_FOR_EACH_WELL_H_
#define BS_FOR_EACH_WELL_H_

namespace blue_sky {

  /**
   * \brief  For each well in facility_list calls method method
   * \param  facility_list
   * \param  method
   * */
  template <typename method_t>
  void
  for_each_well (facility_manager &facility_list, method_t method)
  {
    facility_manager::well_const_iterator_t wb = facility_list.wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list.wells_end ();
    typedef well well_t;
    typedef smart_ptr <well_t, true> sp_well_t;

    for (; wb != we; ++wb)
      {
        sp_well_t well (wb->second, bs_dynamic_cast ());
        BS_ASSERT (well);

        if (!well->is_no_connections ())
          {
            method (*well);
          }
      }
  }

  /**
   * \brief  For each facility in facility_list calls method method
   * \param  facility_list
   * \param  method
   * */
  template <typename method_t>
  void
  for_each_facility (facility_manager &facility_list, method_t method)
  {
    facility_manager::well_const_iterator_t wb = facility_list.wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list.wells_end ();

    for (; wb != we; ++wb)
      {
        method (*wb->second);

        //if (well->get_connections_count ())
        //  {
        //    method (*well);
        //  }
      }
  }

  /**
   * \brief  For each well in facilit_list calls method method 
   *         while method not return true
   * \param  facility_list
   * \param  method
   * \return True if for any well in list method returns true
   * */
  template <typename method_t>
  bool
  for_each_well_while_cond (facility_manager &facility_list, method_t method)
  {
    facility_manager::well_const_iterator_t wb = facility_list.wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list.wells_end ();
    typedef well well_t;
    typedef smart_ptr <well_t, true> sp_well_t;

    for (; wb != we; ++wb)
      {
        sp_well_t well (wb->second, bs_dynamic_cast ());
        BS_ASSERT (well);
        if (!well->is_no_connections () && method (*well))
          {
            return true;
          }
      }

    return false;
  }
  /**
   * \brief  Accumulates values returned by method for each well in 
   *         facility_list
   * \param  facility_list
   * \param  method
   * \return Accumulated value
   * */
  template <typename method_t>
  size_t
  for_each_well_acc (facility_manager &facility_list, method_t method)
  {
    facility_manager::well_const_iterator_t wb = facility_list.wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list.wells_end ();
    typedef well well_t;
    typedef smart_ptr <well_t, true> sp_well_t;

    size_t acc = 0;
    for (; wb != we; ++wb)
      {
        sp_well_t well (wb->second, bs_dynamic_cast ());
        BS_ASSERT (well);

        if (!well->is_no_connections ())
          {
            acc += method (*well);
          }
      }

    return acc;
  }

} // namespace blue_sky

#endif  // #ifndef BS_FOR_EACH_WELL_H_

