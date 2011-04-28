/**
 * \file hdf5_functions.h
 * \brief Functions for easy writing/reading data to/from hdf5 file
 * \author Vadim Alimguzhin
 * \date 09.10.2008
 */

#ifndef HDF5_FUNCTIONS_H_
#define HDF5_FUNCTIONS_H_

#include "hdf5.h"
#include "hdf5_type.h"

namespace blue_sky
  {

  /**
   * \brief Write scalar data to hdf5 file at specified location
   * \param location -- location in hdf5 file where to write
   * \param name -- dataset name
   * \param data -- data to write
   */
  template<typename T> int
  write_data_to_hdf5(const hid_t &location, const std::string &name, const T &data)
  {
    const hid_t &datatype = get_hdf5_type(data);
    hid_t dataspace = H5Screate(H5S_SCALAR);
    hid_t cparms = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(cparms, H5D_COMPACT);
    hid_t dataset = H5Dcreate(location, name.c_str(), datatype, dataspace, cparms);
    if (dataset < 0)
      {
        printf ("Error: can't create dataset %s while writing hdf5!\n", name.c_str());
        return -1;
      }
    hid_t fspace = H5Dget_space(dataset);
    herr_t status = H5Dwrite(dataset, datatype, NULL, fspace, NULL, &data);

    H5Sclose (dataspace);
    H5Sclose (fspace);
    H5Pclose (cparms);
    H5Dclose (dataset);
    return (int) status;
  }

  /**
   * \brief Read scalar data from hdf5 file at specified location
   * \param location -- location in hdf5 file from where to read
   * \param name -- dataset name
   * \param data
   */
  template<typename T> int
  read_data_from_hdf5(const hid_t &location, const std::string &name, T &data)
  {
    const hid_t &datatype = get_hdf5_type(data);
    htri_t status = H5Lexists(location, name.c_str(), NULL);
    if (!status)
      {
        printf ("Error: can't open dataset %s while reading hdf5 (not exist)!\n", name.c_str());
        return -1;
      }
    hid_t dataset = H5Dopen(location, name.c_str());
    herr_t status2 = H5Dread(dataset, datatype, NULL, NULL, H5P_DEFAULT, &data);
    H5Dclose(dataset);
    return (int) status2;
  }

  /**
   * \brief Write md array to hdf5 file at specified location
   * \param location -- location in hdf5 file where to write
   * \param name -- data name
   * \param data -- array to write
   * \param rank -- rank of the array
   * \param dims -- array's dimensions
   */
  template<typename T> int
  write_data_to_hdf5(const hid_t &location, const std::string &name, const T *data, int rank, const hsize_t *dims)
  {
    const hid_t &datatype = get_hdf5_type(data);
    hid_t dataspace = H5Screate_simple(rank, dims, NULL);
    hid_t dataset = H5Dcreate(location, name.c_str(), datatype, dataspace, H5P_DEFAULT);
    if (dataset < 0)
      {
        printf ("Error: can't create dataset %s while writing hdf5!\n", name.c_str());
        return -1;
      }
    herr_t status = H5Dwrite(dataset, datatype, 0, 0, H5P_DEFAULT, data);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return (int) status;
  }

  /**
   * \brief Write 1d array to hdf5 file at specified location
   * \param location -- location in hdf5 file where to write
   * \param name -- data name
   * \param data -- array to write
   * \param length -- array's length
   */
  template<typename T> int
  write_data_to_hdf5(const hid_t &location, const std::string &name, const T *data, int length)
  {
    hsize_t dims[] = {length};
    return write_data_to_hdf5(location, name, data, 1, dims);
  }

  /**
   * \brief Read array from hdf5 file at specified location
   * \brief memory is allocated inside function!
   * \param location -- location in hdf5 file from where to read
   * \param name -- data name
   * \param data
   * \return number of read data elements
   */
  template<typename T> int
  read_data_from_hdf5(const hid_t &location, const std::string &name, T **data)
  {
    const hid_t &datatype = get_hdf5_type<T>();
    htri_t status = H5Lexists(location, name.c_str(), NULL);
    if (!status)
      {
        printf ("Error: can't open dataset %s while reading hdf5 (not exist)!\n", name.c_str());
        return -1;
      }
    hid_t dataset = H5Dopen(location, name.c_str());
    hid_t dataspace  = H5Dget_space(dataset);
    int array_length = (int) H5Sget_simple_extent_npoints(dataspace);
    if (*data)
      delete[] *data;
    *data = new T[array_length];
    herr_t status2 = H5Dread(dataset, datatype, NULL, NULL, H5P_DEFAULT, *data);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return (int) status2;//array_length;
  }

  /**
   * \brief Read array with specified length from hdf5 file at specified location
   * \brief memory is allocated inside function!
   * \param location -- location in hdf5 file from where to read
   * \param name -- data name
   * \param data
   * \param length - array length
   * \return 0 if success
   */
  template<typename T> int
  read_data_from_hdf5(const hid_t &location, const std::string &name, T **data, int length)
  {
    const hid_t &datatype = get_hdf5_type<T>();
    htri_t status = H5Lexists(location, name.c_str(), NULL);
    if (!status)
      {
        printf ("Error: can't open dataset %s while reading hdf5 (not exist)!\n", name.c_str());
        return -1;
      }
    hid_t dataset = H5Dopen(location, name.c_str());
    hid_t dataspace  = H5Dget_space(dataset);
    int array_length = (int) H5Sget_simple_extent_npoints(dataspace);
    if (length != array_length)
      return -1;
    if (*data)
      delete[] *data;
    *data = new T[array_length];
    herr_t status2 = H5Dread(dataset, datatype, NULL, NULL, H5P_DEFAULT, *data);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return (int) status2;
  }

  /**
   * \brief Read array with specified length from hdf5 file at specified location
   * \param location -- location in hdf5 file from where to read
   * \param name -- data name
   * \param data
   * \param length - array length
   * \return 0 if success
   */
  template<typename T> int
  read_data_from_hdf5(const hid_t &location, const std::string &name, T *data, int length)
  {
    const hid_t &datatype = get_hdf5_type<T>();
    htri_t status = H5Lexists(location, name.c_str(), NULL);
    if (!status)
      {
        printf ("Error: can't open dataset %s while reading hdf5 (not exist)!\n", name.c_str());
        return -1;
      }
    hid_t dataset = H5Dopen(location, name.c_str());
    hid_t dataspace  = H5Dget_space(dataset);
    int array_length = (int) H5Sget_simple_extent_npoints(dataspace);
    if (length != array_length)
      return -1;
    herr_t status2 = H5Dread(dataset, datatype, NULL, NULL, H5P_DEFAULT, data);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return (int) status2;
  }

} // namespace blue_sky

#endif // #ifndef HDF5_FUNCTIONS_H_

