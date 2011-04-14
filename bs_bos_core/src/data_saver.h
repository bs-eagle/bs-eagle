/**
 * */
#ifndef BS_EAGLE_DATA_SAVER_H_
#define BS_EAGLE_DATA_SAVER_H_

#include "calc_model.h"
#include "rs_mesh_iface.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "jac_matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "calc_well.h"
#include "facility_manager.h"

#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace blue_sky {

  struct BS_API_PLUGIN data_saver
  {
    typedef t_double          item_t;

    typedef calc_model        calc_model_t;
    typedef rs_mesh_iface     mesh_iface_t;
    typedef jac_matrix_iface  jacobian_matrix_t; 
    typedef facility_manager  facility_manager_t;

    typedef facility_manager_t::well_const_iterator_t      well_iterator_t;

    typedef smart_ptr <calc_model_t>      sp_calc_model_t;
    typedef smart_ptr <mesh_iface_t>      sp_mesh_iface_t;
    typedef smart_ptr <jacobian_matrix_t> sp_jmatrix_t;

  public:
    /**
     * \brief  Constructs data saver
     * */
    data_saver ();

    /**
     * \brief  Opens data storage
     * \param  name
     * */
    void
    open_storage (const std::string &name);

    /**
     * \brief  Writes wells and perforation data to storage
     * \param  calc_model
     * */
    void
    write_well_results (const sp_calc_model_t &calc_model,
      well_iterator_t wb, 
      const well_iterator_t &we, 
      double time);

    /**
     * \brief  Writes FIP data to storage
     * \param  calc_model
     * */
    void
    write_fip_results (const sp_calc_model_t &calc_model);

    /**
     * \brief  Writes calc_model data to storage
     * \param  calc_model
     * \param  jmx Jacobian matrix (to save CFL data)
     * \param  large_time_step_num Number of current large time-step
     * \param  total_time_step_num
     * \param  time
     * */
    void
    write_calc_model_data (const sp_calc_model_t &calc_model,
      const sp_jmatrix_t &jmx,
      size_t large_time_step_num,
      size_t total_time_step_num,
      double time);

    /**
     * \brief  Writes inital status of mesh to storage
     * \param  mesh
     * */
    void
    write_mesh (const sp_mesh_iface_t &mesh);

    /**
     * \brief  Writes starting date to storage
     * \param  date
     * */
    void
    write_starting_date (const boost::posix_time::ptime &date);

  private:

    struct impl;

    boost::shared_ptr <impl>    impl_; //!< data_saver implementation
  };

} // namespace blue_sky


#endif // #ifndef BS_EAGLE_DATA_SAVER_H_

