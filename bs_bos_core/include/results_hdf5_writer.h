/**
 *       \file  results_hdf5_writer.h
 *      \brief  Writes [well|fip]_results to HDF5 storage
 *     \author  Original author Sayfulling Ilshat, splited to
 *              separate file by Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  25.08.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_BOS_CORE_RESULTS_HDF5_WRITER_H_
#define BS_BOS_CORE_RESULTS_HDF5_WRITER_H_

#ifdef _HDF5
#include "bs_hdf5_storage.h"
//#include "well_results_storage.h"
//#include "fip_results_storage.h"
#include "rs_smesh_iface.h"

namespace blue_sky
  {
  namespace hdf5
    {

    ///**
    // * \brief Writes well_results_storage data to hdf5 file
    // * \param hdf5 File to write to
    // * \param wres well_results_storage
    // * \param write_conn_data if true then connections data will be written to file
    // */
    //void write_well_results (const bs_hdf5_storage &hdf5, smart_ptr <well_results_storage> &well_res, bool write_conn_data);

    ///**
    // * \brief Writes fip_results_storage data to hdf5 file
    // * \param hdf5 File to write to
    // * \param fres fip_results_storage
    // */
    //void write_fip_results (const bs_hdf5_storage &hdf5, smart_ptr <fip_results_storage> &fip_res);

    /**
     * \brief Writes mesh data to hdf5 file
     * \param hdf5 File to write to
     * \param mesh
     * */
    void write_mesh_to_hdf5 (const smart_ptr <bs_hdf5_storage, true> &hdf5, const smart_ptr <rs_smesh_iface, true> &mesh);

    namespace detail
      {
      /**
       * \brief  Auxillary function
       * */
      template <typename params_t>
      void add_params (const bs_hdf5_storage &hdf5, const hid_t &hid, const std::string &name, const params_t &params, size_t size);

      /**
       * \brief  Auxillary function
       * */
      template <typename dates_t>
      void add_dates (const bs_hdf5_storage &hdf5, const hid_t &hid, const dates_t &dates);
    } //namespace detail

  } // namespace hdf5
} // namespace blue_sky

#endif // #ifdef _HDF5

#endif // #ifndef BS_BOS_CORE_RESULTS_HDF5_WRITER_H_

