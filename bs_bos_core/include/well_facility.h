/**
 *       \file  well_facility.h
 *      \brief  Base interface for well facilities
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  07.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_EAGLE_WELL_FACILITY_H_
#define BS_EAGLE_WELL_FACILITY_H_

namespace blue_sky {

  class well;

namespace wells {

  /**
   * \class well_facility_iface
   * \brief Base interface for well facilities
   * */
  struct well_facility_iface : public objbase
  {
    virtual ~well_facility_iface () {}

    /**
     * \brief  Performs actions before well process function
     * \param  well
     * */
    virtual void
    pre_process (well *well) = 0;

    /**
     * \brief  Performs actions after well process function
     * \param  well
     * */
    virtual void
    post_process (well *well) = 0;
  };

} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_EAGLE_WELL_FACILITY_H_

