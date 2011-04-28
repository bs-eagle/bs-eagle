#ifndef __BS_HDF5_STORAGE_H
#define __BS_HDF5_STORAGE_H

#include "bs_abstract_storage.h"//TODO: remove this by #include "bs_object_base.h"
#include "hdf5_functions.h"
#include "hdf5.h"
#include "throw_exception.h"

namespace blue_sky
  {

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
        htri_t status = H5Lexists (location, name.c_str (), 0);
        return status > 0;
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

  class BS_API_PLUGIN bs_hdf5_storage : public objbase //public bs_abstract_storage
    {
      BLUE_SKY_TYPE_DECL (bs_hdf5_storage);

    private:
      // hdf5 file descriptor
      hid_t file_id;
      // flag is hdf5 file already opened
      int open_flag;

    public:

      typedef std::vector<double> dates_type;
      typedef std::vector<float> d_params_internal_type;
      typedef std::vector<d_params_internal_type> d_params_type;
      typedef std::vector<int> i_params_internal_type;
      typedef std::vector<i_params_internal_type> i_params_type;

      // COMMON FUNCTIONS:


      // Return file descriptor
      hid_t get_file_id () const
      {
        return file_id;
      };

      // Open file
      int open (const char* filename, int flag = ST_ACC_CREATE);

      // Close file
      int close ();

      // Returns true if opened, false else.
      bool is_opened () const;

      // Create group
      hid_t begin_object (const hid_t location, const std::string &name) const;
      hid_t begin_object (const std::string &name) const;

      // Get rank of data set. Returns number of dimension
      int get_rank (const hid_t location, const std::string &name) const;

      // Get dimensions
      int get_dimensions (const hid_t location, const std::string &name, int *dimensions) const;

      // Get name of object in group by it's index
      std::string get_obj_name_by_idx (const hid_t &location, const hsize_t idx) const;

      // read any data from dataset with path of hdf5 file to buf
      template<typename T>
      int get_data (const char *path, T *buf) const;

      // get size of dataset elements with path from file_in
      int get_datatype_size (const char *path) const;

      // flush buffer associated with hdf5 file to disk
      int flush ();

    public:
      int create_group (hid_t location, const char *group_name) const;

      // writing 1D-array to hdf5 file
      template<typename T>
      int write_array (const char *path, const char *dataset_name, const T *data, int n_elements) const;

      // writing 1D-array (as extendable dataset) to hdf5 file
      template<typename T>
      int add_to_results (std::string path, const T *value, int n_elem_global, int n_elem_local, int mpi_offset, double t) const;

      int write_string_to_hdf5 (const hid_t location, const char *dataset_name, const std::string &str) const;
    };

  // check group exists, create if doesn't exist, open group else
  inline int
  bs_hdf5_storage::create_group (const hid_t location, const char *group_name) const
  {
    hid_t group;
    if (!H5Lexists (location, group_name, 0)) // not exist
    group = H5Gcreate (location, group_name, 0);
    else
      group = H5Gopen (location, group_name);
      if (group < 0)
        {
          throw bs_exception ("bs_hdf5_storage", "group_creating_error");
          }
    return group;
  }

  // calls H5Gget_objname_by_idx and return as std::string
  inline std::string
  bs_hdf5_storage::get_obj_name_by_idx (const hid_t &location, const hsize_t idx) const
  {
    char name[100];
    int name_size = H5Gget_objname_by_idx (location, idx, NULL, 0) + 1;//get size of the name, +1 - hdf5 specification
    H5Gget_objname_by_idx (location, idx, name, name_size);//get name
    return std::string (name);
  }

  // writing 1D-array of type T to HDF5 file
  template<typename T> int
  bs_hdf5_storage::write_array (const char *path, const char *dataset_name,
                                const T *data, int n_elements) const
  {
    const int RANK = 1;
    // open/create group
    hid_t mesh_group = create_group (file_id, path);
    // Create dataset
    hsize_t dim[] =
      {
        n_elements
      };
    hid_t space = H5Screate_simple (RANK, dim, NULL);
    hid_t dataset = H5Dcreate (mesh_group, dataset_name, get_hdf5_type <T> (), space, H5P_DEFAULT);
    H5Sclose (space);
    // and write it into the file.
    hid_t plist = H5Pcreate (H5P_DATASET_XFER);
    herr_t status = H5Dwrite (dataset, get_hdf5_type <T> (), NULL, NULL, plist, data);
    status;
    H5Pclose (plist);
    H5Dclose (dataset);
    H5Gclose (mesh_group);
    return 0;
  }

  //! write data of type T and current simulation time
  //! to hdf5 file as extendable dataset
  // std::string path - path in hdf5 file
  // const T *data - data
  // int n_elem_gl - global number of elements (in MPI case)
  // int n_elem_loc - local number of elements (in MPI case)
  // int mpi_offset - start index (in MPI case)
  // double t - current simulation time
  template<typename T> int
  bs_hdf5_storage::add_to_results (std::string path, const T *data, int n_elem_gl, int n_elem_loc, int mpi_offset, double t) const
  {
    hid_t h5_type = get_hdf5_type<T>();
    const int RANK_VALUE = 2;
    const int RANK_TIME = 2;
    const int n_elem_time = 3; // 1. time, 2. max_value, 3. min_value

    // calculate max and min value
    double min = data[0];
    double max = data[0];
    for (int i = 1; i < n_elem_loc; i++)
    {
      min = MIN (min, data[i]);
        max = MAX (max, data[i]);
      }

#ifdef _MPI
    double mpi_data[4]; //2 for send, 2 for recv
    mpi_data[0] = -min; // because max(-a[i]) = -min(a[i])
    mpi_data[1] = max;
    MPI_Reduce (mpi_data, mpi_data + 2, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    min = -mpi_data[3];
    max = mpi_data[4];
#endif //_MPI

    // fill time_data = {time, max, min}
    double time_data[n_elem_time];
    time_data[0] = t;
    time_data[1] = max;
    time_data[2] = min;

    hid_t group, dataset_value, dataset_time, fspace_value, fspace_time, plist_value, plist_time;
    hsize_t size[2];

    hsize_t maxdim_value[2] = {n_elem_gl, H5S_UNLIMITED};
    hsize_t maxdim_time[2] = {n_elem_time, H5S_UNLIMITED};
    hsize_t dims_value_gl[2]  = {n_elem_gl, 1};// dataset dimensions at creation
    hsize_t dims_value_loc[2]  = {n_elem_loc, 1};// dataset dimensions at creation
    hsize_t dims_time[2] = {n_elem_time, 1};
    herr_t status;

    hid_t mspace_time = H5Screate_simple (RANK_TIME, dims_time, maxdim_time);
    hid_t mspace_value = H5Screate_simple (RANK_VALUE, dims_value_loc, maxdim_value);

    std::string time = "dates";
    std::string values = "values";

    if (!H5Lexists (file_id, path.c_str (), NULL)) // if group not exists, create it
    {
      const hsize_t chunk_dims_value[2] =
          {
            n_elem_gl, 1
          };
        const hsize_t chunk_dims_time[2] =
          {
            n_elem_time, 1
          };

        group = H5Gcreate (file_id, path.c_str (), H5P_DEFAULT);

        plist_value = H5Pcreate (H5P_DATASET_CREATE);
        plist_time = H5Pcreate (H5P_DATASET_CREATE);

        H5Pset_chunk (plist_value, RANK_VALUE, chunk_dims_value);
        H5Pset_chunk (plist_time, RANK_TIME, chunk_dims_time);

        fspace_value = H5Screate_simple (RANK_VALUE, dims_value_gl, maxdim_value);
        fspace_time = H5Screate_simple (RANK_TIME, dims_time, maxdim_time);

        dataset_value = H5Dcreate (group, values.c_str (), h5_type, fspace_value, plist_value);
        dataset_time = H5Dcreate (group, time.c_str (), H5T_NATIVE_DOUBLE, fspace_time, plist_time);

        if (dataset_value < 0 || dataset_time < 0)
          {
            throw bs_exception("bs_hdf5_storage", "dataset_create_error in add_to_results!");
            return -2;
          }
        H5Pclose (plist_value);
        H5Pclose (plist_time);
        H5Sclose (fspace_value);
        H5Sclose (fspace_time);

        plist_value = H5Pcreate (H5P_DATASET_XFER);
        plist_time = H5Pcreate (H5P_DATASET_XFER);

        fspace_value = H5Dget_space(dataset_value);
        fspace_time = H5Dget_space(dataset_time);
#ifdef _MPI
        hsize_t offset[2];
        offset[1] = 0;
        offset[0] = mpi_offset; // displasement in arrays with dimension of n_active_elements
        H5Pset_dxpl_mpio (plist_time, H5FD_MPIO_COLLECTIVE);// only one proc
        H5Pset_dxpl_mpio (plist_value, H5FD_MPIO_COLLECTIVE);
        H5Sselect_hyperslab (fspace_value, H5S_SELECT_SET, offset, NULL, dims_value_loc, NULL);
#endif
        status = H5Dwrite (dataset_value, h5_type, mspace_value, fspace_value, plist_value, data);
        status = H5Dwrite (dataset_time, H5T_NATIVE_DOUBLE, mspace_time, fspace_time, plist_time, time_data);
      }
    else // group already exists
      {
        group = detail::open_group (file_id, path.c_str ());
        dataset_value = detail::open_dataset (group, values.c_str ());
        dataset_time = detail::open_dataset (group, time.c_str ());

        hsize_t size_value[2];
        hsize_t size_time[2];

        hid_t space_value  = H5Dget_space (dataset_value);
        hid_t space_time  = H5Dget_space (dataset_time);
        H5Sget_simple_extent_dims (space_value, size_value, NULL);
        H5Sget_simple_extent_dims (space_time, size_time, NULL);
        H5Sclose (space_value);
        H5Sclose (space_time);

        size[0]   = size_value[0];
        size[1]   = size_value[1] + 1;
        H5Dextend (dataset_value, size);

        size[0] = size_time[0];
        size[1] = size_time[1] + 1;
        H5Dextend (dataset_time, size);

        hsize_t     offset[2];
        offset[0] = 0;
        offset[1] = size_value[1];

        plist_value = H5Pcreate (H5P_DATASET_XFER);
#ifdef _MPI
        offset[0] = mpi_offset; // displasement in arrays with dimension of n_active_elements
        H5Pset_dxpl_mpio (plist_value, H5FD_MPIO_COLLECTIVE);
#endif
        fspace_value  = H5Dget_space(dataset_value);
        H5Sselect_hyperslab (fspace_value, H5S_SELECT_SET, offset, NULL, dims_value_loc, NULL);
#ifndef _MPI //there is an error occur, all params seems ok, maybe error in hdf5 lib
        status = H5Dwrite (dataset_value, h5_type, mspace_value, fspace_value, plist_value, data);
#endif
        offset[0] = 0;
        offset[1] = size_time[1];
        plist_time = H5Pcreate (H5P_DATASET_XFER);
        fspace_time  = H5Dget_space (dataset_time);
        H5Sselect_hyperslab (fspace_time, H5S_SELECT_SET, offset, NULL, dims_time, NULL);
        status = H5Dwrite (dataset_time, H5T_NATIVE_DOUBLE, mspace_time, fspace_time, plist_time, time_data);
      }

    H5Sclose (mspace_value);
    H5Sclose (mspace_time);
    H5Sclose (fspace_value);
    H5Sclose (fspace_time);
    H5Pclose (plist_value);
    H5Pclose (plist_time);
    H5Gclose (group);
    H5Dclose (dataset_value);
    H5Dclose (dataset_time);

    return 0;
  }

  // get type of dataset with path
  inline
  int bs_hdf5_storage::get_datatype_size (const char *path) const
  {
    hid_t dataset = detail::open_dataset (file_id, path);
    if (dataset < 0)
    return -1;
    hid_t datatype = H5Dget_type(dataset);
    size_t data_type_size = H5Tget_size(datatype);
    H5Tclose (datatype);
    H5Dclose (dataset);
    return (int) data_type_size;
  }

  // read any data from dataset with path of hdf5 file to buf
  template<typename T> int
  bs_hdf5_storage::get_data (const char *path, T *buf) const
  {
    hid_t h5_type = get_hdf5_type<T>();
    hid_t dataset = detail::open_dataset (file_id, path);
    if (dataset < 0)
    return -2;
    herr_t status = H5Dread(dataset, h5_type, NULL, NULL, H5P_DEFAULT, buf);
    H5Dclose(dataset);
    if (status < 0)
      return -3;
      return 0;
    }

    // AUIXILARY FUNCTIONS

    /*
     *  copy from array of double type to float array, use striding
     *  assume double_data dimension is n_elements * stride
     *  assume float_data dimension is n_elements
     *  for example to copy {y1,y2,y3} from {x1, y1, z1, x2, y2, z2, x2, y2, z2}
     *  set start=1, stride=3
     */
  template<typename T, typename D>
  void copy_to (const T *data, D *float_data, int n_elements, int stride = 1, int start = 0)
  {
    int i, j;
    if (stride >= 1)
      for (i = start, j = 0; i < n_elements; i += stride)
        {
          float_data[j++] = static_cast <D> (data[i]);
        }
    else
      {
        bs_throw_exception (boost::format ("Stride should be greater than 1 (%d)") % stride);
      }
  }

} // namespace blue_sky

#endif
