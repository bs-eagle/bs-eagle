/**
 *       \file  data_saver.cpp
 *      \brief  Saves data to storage
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "data_saver.h"

#ifdef _HDF5
#include "data_saver_hdf5.h"
#else
#include "data_saver_dummy.h"
#endif

#include "rs_smesh_iface.h"
#include "well_connection.h"

namespace blue_sky {

  data_saver::data_saver ()
  {
  }

  void
  data_saver::open_storage (const std::string &name)
  {
    impl_.reset (new impl (name));
  }

  void
  data_saver::write_well_results (const sp_calc_model_t &calc_model, 
    well_iterator_t wb, 
    const well_iterator_t &we, 
    double time)
  {
    BS_ASSERT (impl_);
    impl_->write_well_results (calc_model, wb, we, time);
  }

  void
  data_saver::write_fip_results (const sp_calc_model_t & /*calc_model*/)
  {
    BS_ASSERT (impl_);
    //! \todo Should be implemented
  }

  void
  data_saver::write_calc_model_data (const sp_calc_model_t &calc_model,
    const BS_SP (jacobian) &jmx,
    size_t large_time_step_num,
    size_t total_time_step_num,
    double time)
  {
    BS_ASSERT (impl_);
    impl_->write_calc_model_data (calc_model, jmx, large_time_step_num, total_time_step_num, time);
  }

  void
  data_saver::write_mesh (const sp_mesh_iface_t &mesh)
  {
    BS_ASSERT (impl_);
    impl_->write_mesh (mesh);
  }

  void
  data_saver::write_starting_date (const boost::posix_time::ptime &date)
  {
    BS_ASSERT (impl_);
    impl_->write_starting_date (date);
  }
} // namespace blue_sky

