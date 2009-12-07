/**
 * */
#ifndef BS_EAGLE_WELL_FACILITY_H_
#define BS_EAGLE_WELL_FACILITY_H_

namespace blue_sky {

  template <typename strategy_t>
  class well;

namespace wells {

  template <typename strategy_t>
  struct well_facility_iface : public objbase
  {
    virtual ~well_facility_iface () {}

    /**
     * \brief  Performs actions before well process function
     * \param  well
     * */
    virtual void
    pre_process (well <strategy_t> *well) = 0;

    /**
     * \brief  Performs actions after well process function
     * \param  well
     * */
    virtual void
    post_process (well <strategy_t> *well) = 0;
  };

} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_EAGLE_WELL_FACILITY_H_

