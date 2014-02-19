#ifndef HDF5_FUNCTIONS_H_72b44b6c_85e1_11e0_a55f_2f89c88b9605
#define HDF5_FUNCTIONS_H_72b44b6c_85e1_11e0_a55f_2f89c88b9605
/**
 *       \file  hdf5_functions.h
 *      \brief  Wrappers over HDF5 native function
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *    \comment  No code by Vadim Alimguzhin left, so I remove refs to him
 * */

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "hdf5.h"
#include "hdf5_type.h"
#include "hdf5_type_to_hid.hpp"

#include "throw_exception.h"

#include <boost/tokenizer.hpp>


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
    const hid_t &datatype = get_hdf5_type <hdf5_type_helper <T>::type> ();
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
    const hid_t &datatype = get_hdf5_type <hdf5_type_helper <T>::type> ();
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
    const hid_t &datatype = get_hdf5_type <hdf5_type_helper <T>::type> ();
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
    hsize_t dims[] = {hsize_t(length)};
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
    const hid_t &datatype = get_hdf5_type <hdf5_type_helper <T>::type> ();
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
    const hid_t &datatype = get_hdf5_type <hdf5_type_helper <T>::type> ();
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
    const hid_t &datatype = get_hdf5_type <hdf5_type_helper <T>::type> ();
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

  namespace detail {

    /**
     * \brief  Checks is object name is exists in location
     * \param  location HDF5 object id
     * \param  name Name of the object to check
     * \return True if object exists
     * */
    inline bool 
    is_object_exists (const hid_t &location, const std::string &name) 
    {
      typedef boost::tokenizer <boost::char_separator <char> > tokenizer_t;
      boost::char_separator <char> s ("/");
      tokenizer_t t (name, s);
      std::string path = "";
      for (tokenizer_t::iterator it = t.begin (), e = t.end (); it != e; ++it)
        {
          path += *it;
          htri_t status = H5Lexists (location, path.c_str (), 0);
          if (status <= 0)
            return false;

          path += "/";
        }

      return true;
    }

    /** 
     * \brief Checks if a file with given name is HDF5 file
     * \param name name of file
     * \return True if file if HDF5 file
     * */
    inline bool
    is_hdf5_file (std::string const &name)
    {
      htri_t r = H5Fis_hdf5 (name.c_str ());
      return r > 0;
    }

    /**
     * \brief  Returns object name
     * \param  location HDF5 object id
     * \return Name of the object on success otherwise throws exception
     * */
    inline std::string
    object_name (hid_t location, const std::string &path)
    {
      char *name = 0;
      int name_length = (int) H5Iget_name (location, name, 0);
      if (name_length < 0)
        {
          bs_throw_exception (boost::format ("Can't get HDF5 object name, path: %s") % path);
        }

      name = new char[name_length + 1];
      memset (name, 0, name_length + 1);
      H5Iget_name (location, name, name_length + 1);

      std::string n (name, name_length);
      delete [] name;

      return n;
    }

    /**
     * \brief  Open group if exists, if no throws exception
     * \param  location HDF5 object id
     * \param  path Path of the group to open
     * \return Group id on success otherwise throws exception
     * */
    inline hid_t
    open_group (const hid_t location, const char *path) 
    {
      if (!H5Lexists (location, path, 0)) // not exist
        {
          bs_throw_exception (boost::format ("Group %s/%s is not exists") % object_name (location, path) % path);
        }

      hid_t group = H5Gopen (location, path);
      if (group < 0)
        {
          bs_throw_exception (boost::format ("Can't open group %s/%s") % object_name (location, path) % path);
        }

      return group;
    }

    /**
     * \brief  Open dataset if exists, if no throws exception
     * \param  location HDF5 object id
     * \param  path Path of the dataset to open
     * \return Dataset id on success otherwise throws exception
     * */
    inline hid_t
    open_dataset (const hid_t location, const char *path) 
    {
      if (!H5Lexists(location, path, 0)) // not exist
        {
          bs_throw_exception (boost::format ("Dataset %s/%s is not exists") % object_name (location, path) % path);
        }

      hid_t ds = H5Dopen (location, path);
      if (ds < 0)
        {
          bs_throw_exception (boost::format ("Can't open dataset %s/%s") % object_name (location, path) % path);
        }

      return ds;
    }
  }


} // namespace blue_sky

#endif // #ifndef HDF5_FUNCTIONS_H_

