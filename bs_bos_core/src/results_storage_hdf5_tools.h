/**
 * \file results_storage_hdf5_tools.h
 * \brief
 * \author Sergey Miryanov
 * \date 25.08.2009
 * */
#ifndef BS_BOS_CORE_RESULTS_STORAGE_HDF5_TOOLS_H_
#define BS_BOS_CORE_RESULTS_STORAGE_HDF5_TOOLS_H_

#ifdef _HDF5

namespace blue_sky {
namespace hdf5 {
namespace detail {

  template <typename dates_t>
  inline void
  add_dates (bs_hdf5_storage &hdf5, const hid_t &hid, const dates_t &dates)
  {
    hid_t dataset;
    if (!hdf5.is_object_exists (hid, "dates"))
      {
        // creating dataset_dates
        hsize_t dims[]      = {0};
        hsize_t dims_max[]  = {H5S_UNLIMITED};

        hid_t plist         = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dataspace     = H5Screate_simple(1, dims, dims_max);

        // set the dataset to be chunked
        hsize_t chunk_dims  = 1;
        H5Pset_chunk(plist, 1, &chunk_dims);

        dataset = H5Dcreate(hid, "dates", H5T_NATIVE_DOUBLE, dataspace, plist);
        H5Pclose(plist);
        H5Sclose(dataspace);
      }
    else
      {
        dataset = hdf5.open_dataset (hid, "dates");
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
    htri_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mspace, fspace, plist, &dates[0]);
    H5Pclose(plist);

    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(dataset);
  }

  template <typename params_t>
  inline void
  add_params (bs_hdf5_storage &hdf5, const hid_t &hid, const std::string &name, const params_t &params, size_t size)
  {
    hid_t dataset;
    if (!hdf5.is_object_exists (hid, name))
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
        dataset = hdf5.open_dataset (hid, name.c_str ());
      }

    // determine new dims of dataset_d_params and extend it
    hsize_t dims_old[2]   = {0};
    hsize_t dims_memory[] = {params.size()};

    hid_t fspace  = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(fspace, dims_old, NULL);
    H5Sclose(fspace);

    hsize_t dims_new[] =
      {
        params.size()
        , dims_old[1] + size
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
      }

    H5Pclose(plist);
    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(dataset);
  }

  inline void
  create_group (bs_hdf5_storage &hdf5, const hid_t &hid, const std::string &name, const std::string &sub_name)
  {
    // write data for the current well
    hid_t dataspace_group = H5Screate(H5S_SCALAR);
    hid_t datatype_group  = H5Tcopy (H5T_C_S1); // create string of given length
    htri_t status         = H5Tset_size (datatype_group, sub_name.length());
    if (status < 0)
      {
        bs_throw_exception (boost::format ("Can't set size for datatype of group, sub_name [%s]") % sub_name);
      }
    hid_t dataset_group   = H5Dcreate(hid, name.c_str (), datatype_group, dataspace_group, H5P_DEFAULT);
    if (dataset_group >= 0)//succesfully created
      {
        hid_t fspace  = H5Dget_space(dataset_group);
        hid_t plist   = H5Pcreate(H5P_DATASET_XFER);
        status        = H5Dwrite(dataset_group, datatype_group, dataspace_group, fspace, plist, &sub_name[0]);

        H5Pclose(plist);
        H5Sclose(fspace);
        H5Dclose(dataset_group);
      }
    else
      {
        bs_throw_exception ("Can't create dataset_group");
      }

    H5Tclose(datatype_group);
  }

} // namespace detail
} // namespace hdf5
} // namespace blue_sky

#endif // #ifdef _HDF5

#endif  // #ifdef BS_BOS_CORE_RESULTS_STORAGE_HDF5_TOOLS_H_

