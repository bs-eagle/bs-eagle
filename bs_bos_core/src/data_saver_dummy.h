/**
 * */
#ifndef BS_EAGLE_DATA_SAVER_DUMMY_H_
#define BS_EAGLE_DATA_SAVER_DUMMY_H_

namespace blue_sky {

  struct data_saver::impl
  {
    typedef t_double   item_t;

    typedef calc_model        calc_model_t;
    typedef rs_mesh_iface     mesh_iface_t;
    typedef facility_manager  facility_manager_t;

    typedef facility_manager_t::well_const_iterator_t      well_iterator_t;

    typedef smart_ptr <calc_model_t>      sp_calc_model_t;
    typedef smart_ptr <mesh_iface_t>      sp_mesh_iface_t;

    impl (const std::string &/*name*/)
    {
      BOSOUT (section::save_data, level::low)
        << "Default implementation of data_saver used"
        << bs_end;
    }

    void
    write_well_results (const sp_calc_model_t &/*calc_model*/,
      well_iterator_t /*wb*/,
      const well_iterator_t &/*we*/,
      double /*time*/)
    {
    }

    void
    write_calc_model_data (const sp_calc_model_t &/*calc_model*/,
      const BS_SP (jacobian) &/*jmx*/,
      size_t /*large_time_step_num*/,
      size_t /*total_time_step_num*/,
      double /*time*/)
    {
    }

    void
    write_mesh (const sp_mesh_iface_t &/*mesh*/)
    {
    }

    void
    write_starting_date (const boost::posix_time::ptime &/*date*/)
    {
    }
  };

} // namespace blue_sky

#endif // #ifndef BS_EAGLE_DATA_SAVER_DUMMY_H_

