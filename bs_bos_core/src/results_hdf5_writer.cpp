/**
 *       \file  results_hdf5_writer.cpp
 *      \brief  Writes [well|fip]_results to HDF5 storage
 *     \author  Original author Sayfulling Ilshat, splited to
 *              separate file by Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  25.08.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "results_hdf5_writer.h"

#ifdef _HDF5

namespace blue_sky
  {
  namespace hdf5
    {

    void
    write_fip_results (const bs_hdf5_storage &hdf5, smart_ptr <fip_results_storage> &fres)
    {
      if (!fres)
        {
          BOSWARN (section::save_data, level::warning) << "Fip results storage is null" << bs_end;
          return ;
        }

      hid_t group_fip = hdf5.begin_object ("fip");

      fip_results_storage::fip_type::iterator iter, e;
      for (iter = fres->fip.begin (), e = fres->fip.end (); iter != e; ++iter)
        {
          char group_name[100];
          sprintf (group_name, "%d", (*iter).first);

          // group for fip: check if exist, then open, create else
          hid_t group = hdf5.begin_object (group_fip, group_name);

          // write data for the current fip_results
          int dates_size = (int) iter->second.dates.size();
          if (dates_size > 0)
            {
              detail::add_dates (hdf5, group, iter->second.dates);
              detail::add_params (hdf5, group, "d_params", iter->second.d_params, dates_size);
            }
          H5Gclose (group);
        }

      H5Gclose (group_fip);
    }

    void
    write_well_results (const bs_hdf5_storage &hdf5, smart_ptr <well_results_storage> &wres, bool write_conn_data)
    {
      if (!wres)
        {
          BOSWARN (section::save_data, level::warning) << "Well results storage is null" << bs_end;
          return ;
        }

      hid_t group_wells, group_well;
      hid_t group_connections, group_connection;

#ifdef _MPI
      //MPI: each process writes wells/conns to own file
      file = &indiv_file;
      group_wells = *file;
#else
      // check if exist, then open, create else
      group_wells = hdf5.begin_object ("wells");
#endif
      well_results_storage::wells_type::iterator iter, e;
      for (iter = wres->wells.begin(), e = wres->wells.end(); iter != e; ++iter)
        {
#ifdef _MPI
          // in mpi decomposition we operate with global mesh in each proc
          // therefore all wells are present in each proc
          // but because other processor's domains are inactive
          // connections of wells, which are in those domains not exists
          // for writing only own wells we use [well.connections.size() > 0] condition
          // but own wells with no connections will be skipped too :(
          // we also send param i_need_write to methods for writing connections,
          // but it's not neccessary, because loop will be empty ;)
          if (!iter->second.connections.size())
            {
              continue;
            }
#endif

          group_well = hdf5.begin_object (group_wells, iter->first);

          // write well's group name
          hdf5.write_string_to_hdf5 (group_well, "group", iter->second.group);

          int dates_size = (int) iter->second.dates.size();
          if (dates_size > 0)
            {
              detail::add_dates (hdf5, group_well, iter->second.dates);
              detail::add_params (hdf5, group_well, "d_params", iter->second.d_params, dates_size);
              detail::add_params (hdf5, group_well, "i_params", iter->second.i_params, dates_size);

              iter->second.clear ();
            }

          // write connections data
          if (write_conn_data)
            {
              // group for connection: check if exist, then open, create else
              group_connections = hdf5.begin_object (group_well, "connections");

              // write connections
              char connection_name[100];
              well_results::conn_type::iterator iter_conn;
              for (iter_conn = iter->second.connections.begin(); iter_conn != iter->second.connections.end(); ++iter_conn)
                {
                  // create group for connection
                  sprintf (connection_name, "%d", iter_conn->first);

                  // group for connection: check if exist, then open, create else
                  group_connection = hdf5.begin_object (group_connections, connection_name);

                  int cdates_size = (int) iter_conn->second.dates.size();
                  if (cdates_size > 0)
                    {
                      detail::add_dates (hdf5, group_connection, iter_conn->second.dates);
                      detail::add_params (hdf5, group_connection, "d_params", iter_conn->second.d_params, cdates_size);
                      detail::add_params (hdf5, group_connection, "i_params", iter_conn->second.i_params, cdates_size);

                      iter_conn->second.clear ();
                    }
                  H5Gclose (group_connection);
                } // conn loop
            } // if (write_conn_data)
          H5Gclose (group_well);
        }//wells loop
#ifndef _MPI
      H5Gclose (group_wells);
#endif
    }

    template<class strategy_t>
    void write_mesh_to_hdf5 (const smart_ptr <bs_hdf5_storage, true> &hdf5, const smart_ptr <rs_smesh_iface<strategy_t>, true> &mesh)
    {
      const int INITIAL_DATA_SIZE = 5;
      int values[INITIAL_DATA_SIZE];

      typename rs_smesh_iface<strategy_t>::index_point3d_t dims = mesh->get_dimens ();

      values[0] = dims[0];
      values[1] = dims[1];
      values[2] = dims[2];
      values[3] = mesh->get_n_active_elements ();
      values[4] = mesh->get_n_active_elements ();//mesh->n_planes

      hdf5->write_array ("/mesh","initial_data", values, 5);
      hdf5->write_array ("/mesh","original_elements", &(mesh->get_ext_to_int())[0], (int)mesh->get_ext_to_int().size());
      hdf5->write_array ("/mesh","original_planes", &(mesh->get_ext_to_int())[0], (int)mesh->get_ext_to_int().size());//TODO
    }

    namespace detail
      {

      template <typename dates_t>
      void
      add_dates (const bs_hdf5_storage &hdf5, const hid_t &hid, const dates_t &dates)
      {
        hid_t dataset;
        if (!blue_sky::detail::is_object_exists (hid, "dates"))
          {
            // creating dataset_dates
            hsize_t dims[]      = {0};
            hsize_t dims_max[]  = {H5S_UNLIMITED};

            hid_t plist         = H5Pcreate(H5P_DATASET_CREATE);
            hid_t dataspace     = H5Screate_simple(1, dims, dims_max);

            // set the dataset to be chunked
            hsize_t chunk_dims  = 1;
            H5Pset_chunk(plist, 1, &chunk_dims);

            dataset = H5Dcreate(hid, "dates", get_hdf5_type (dates), dataspace, plist);
            H5Pclose(plist);
            H5Sclose(dataspace);
          }
        else
          {
            dataset = blue_sky::detail::open_dataset (hid, "dates");
          }

        // determine new dims of dataset_dates and extend it
        hsize_t dims_old;
        hsize_t dims_memory = dates.size();

        hid_t fspace  = H5Dget_space(dataset);
        H5Sget_simple_extent_dims(fspace, &dims_old, NULL);
        hsize_t dims_new = dims_old + dims_memory;
        H5Dextend(dataset, &dims_new);

        // creating dataspaces
        hid_t mspace = H5Screate_simple(1, &dims_memory, NULL);
        fspace = H5Dget_space(dataset);

        hsize_t count = dims_memory;
        hsize_t start = dims_old;
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &start, NULL, &count, NULL);

        // writing
        hid_t plist = H5Pcreate(H5P_DATASET_XFER);
        htri_t status = H5Dwrite(dataset, get_hdf5_type (dates), mspace, fspace, plist, &dates[0]);
        status;
        H5Pclose(plist);

        H5Sclose(mspace);
        H5Sclose(fspace);
        H5Dclose(dataset);
      }

      template <typename params_t>
      void
      add_params (const bs_hdf5_storage &hdf5, const hid_t &hid, const std::string &name, const params_t &params, size_t size)
      {
        hid_t dataset;
        if (!blue_sky::detail::is_object_exists (hid, name))
          {
            // creating dataset_d_params
            hsize_t dims[]      = {params.size(), 0};
            hsize_t dims_max[]  = {params.size(), H5S_UNLIMITED};
            hid_t dataspace     = H5Screate_simple(2, dims, dims_max);

            // set the dataset to be chunked
            hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk_dims[] = {params.size(), 1};
            H5Pset_chunk(plist, 2, chunk_dims);
            dataset = H5Dcreate(hid, name.c_str (), get_hdf5_type <typename params_t::value_type::value_type> (), dataspace, plist);
            H5Pclose(plist);
            H5Sclose(dataspace);
          }
        else
          {
            dataset = blue_sky::detail::open_dataset (hid, name.c_str ());
          }

        // determine new dims of dataset_d_params and extend it
        hsize_t dims_old[2]   =
          {
            0
          };
        hsize_t dims_memory[] =
          {
            params.size()
          };

        hid_t fspace  = H5Dget_space(dataset);
        H5Sget_simple_extent_dims(fspace, dims_old, NULL);
        H5Sclose(fspace);

        hsize_t dims_new[] =
          {
            params.size(), dims_old[1] + size
          };

        H5Dextend(dataset, dims_new);

        // creating dataspaces
        hid_t mspace = H5Screate_simple(1, dims_memory, NULL);
        fspace = H5Dget_space(dataset);
        hsize_t count[] =
          {
            params.size(), 1
          };

        // writing
        hid_t plist = H5Pcreate(H5P_DATASET_XFER);
        typename params_t::value_type buffer;
        buffer.resize(params.size());

        for (size_t i = 0; i < size; i++)
          {
            for (size_t j = 0; j < params.size(); j++)
              {
                typename params_t::value_type::value_type b = 0;
                if (params[j].size () > i)
                  b = params[j][i];

                buffer[j] = b;
              }

            hsize_t start[2] = {0, dims_old[1] + i};
            H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
            htri_t status = H5Dwrite(dataset, get_hdf5_type <typename params_t::value_type::value_type> (), mspace, fspace, plist, &buffer[0]);
            status;
          }

        H5Pclose(plist);
        H5Sclose(mspace);
        H5Sclose(fspace);
        H5Dclose(dataset);
      }

    } // namespace detail

    template void write_mesh_to_hdf5 <base_strategy_fi> (const smart_ptr <bs_hdf5_storage, true> &hdf5, const smart_ptr <rs_smesh_iface<base_strategy_fi>, true> &mesh);
    template void write_mesh_to_hdf5 <base_strategy_di> (const smart_ptr <bs_hdf5_storage, true> &hdf5, const smart_ptr <rs_smesh_iface<base_strategy_di>, true> &mesh);
    template void write_mesh_to_hdf5 <base_strategy_mixi> (const smart_ptr <bs_hdf5_storage, true> &hdf5, const smart_ptr <rs_smesh_iface<base_strategy_mixi>, true> &mesh);

  } // namespace hdf5
} // namespace blue_sky

#endif //#ifdef _HDF5

