#ifndef __BS_HDF5_STORAGE_H
#define __BS_HDF5_STORAGE_H

#include "bs_abstract_storage.h"//TODO: remove this by #include "bs_object_base.h"
#include "hdf5_functions.h"
#include "hdf5.h"
#include "throw_exception.h"

#include "shared_vector.h"
#include "seq_vector.h"
#include "constants.h"

#include <boost/type_traits.hpp>
#include <boost/array.hpp>
#include "auto_value.h"

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
        htri_t status = H5Lexists (location, name.c_str (), NULL);
        return status > 0;
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
        int name_length = (int) H5Iget_name (location, name, NULL);
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
        if (!H5Lexists (location, path, NULL)) // not exist
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
        if (!H5Lexists(location, path, NULL)) // not exist
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
    if (!H5Lexists (location, group_name, NULL)) // not exist
    group = H5Gcreate (location, group_name, NULL);
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

  namespace hdf5 {
  namespace private_ {

    template <typename T>
    const hid_t &get_hdf5_type_ (const T &)
    {
      return get_hdf5_type <T> ();
    }
    template <>
    inline const hid_t &
    get_hdf5_type_ <main_var_type> (const main_var_type &)
    {
      return H5T_NATIVE_UINT;
    }

    template <typename T>
    inline const hid_t &
    get_hdf5_type_ (const seq_vector <T> &)
    {
      return get_hdf5_type_ <T> (T ());
    }


    struct hdf5_buffer__ 
    {
      typedef size_t size_type;

      template <typename T, typename Y>
      hdf5_buffer__ (const Y &, const T *data, size_type size, size_type stride, size_type offset)
      : type (get_hdf5_type_ (Y ()))
      , size (size)
      , owner_ (stride || !boost::is_same <typename boost::remove_cv <T>::type, typename boost::remove_cv <Y>::type> ())
      , data_ (!owner_ ? (void *)data : (void *)(new Y[stride ? (size / stride) : size]))
      {
        if (owner_)
          {
            if (!stride)
              stride = 1;

            size = size / stride;
            for (size_type i = 0, j = 0; i < size; i += stride, ++j)
              {
                static_cast <Y *> (data_)[j] = static_cast <Y> (data[i]);
              }
          }
      }

      hdf5_buffer__ (void *data, size_type size, hid_t type)
      : type (type)
      , size (size)
      , owner_ (false)
      , data_ (data)
      {
      }

      ~hdf5_buffer__ ()
      {
        if (owner_)
          delete [] data_;
      }

      void *
      data () const
      {
        return data_;
      }

    public:
      hid_t       type;
      size_type   size;

      bool        owner_;
      void        *data_;
    };

    template <typename T>
    struct get_size
    {
      enum { static_size = T::static_size, };
    };

    template <>
    struct get_size <size_t>
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <int>
    {
      enum { static_size = 1,  };
    };

    template <>
    struct get_size <float>
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <double>
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <auto_value <double> >
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <auto_value <float> >
    {
      enum { static_size = 1, };
    };

    template <typename T>
    struct get_type
    {
      typedef typename get_type <typename T::lhs_type>::type type;
    };

    template <>
    struct get_type <size_t>
    {
      typedef size_t type;
    };
    template <>
    struct get_type <int>
    {
      typedef int type;
    };
    template <>
    struct get_type <float>
    {
      typedef float type;
    };
    template <>
    struct get_type <double>
    {
      typedef double type;
    };
    template <>
    struct get_type <auto_value <double> >
    {
      typedef double type;
    };
    template <>
    struct get_type <auto_value <float> >
    {
      typedef float type;
    };

    template <size_t s>
    struct buffer_filler 
    {
      template <typename R, typename B>
      static void
      fill (const R &h, B &buffer)
      {
        buffer_filler <R::lhs_type::static_size>::fill (h.lhs, buffer);
        buffer[R::static_size - 1] = h.rhs;
      }
    };
    template <>
    struct buffer_filler <2>
    {
      template <typename R, typename B>
      static void
      fill (const R &h, B &buffer)
      {
        buffer[0] = h.lhs;
        buffer[1] = h.rhs;
      }
    };

    struct hdf5_pod__
    {
      typedef size_t size_type;

      hdf5_pod__ ()
      : type (0)
      , size (0)
      , data_ (0)
      {
      }

      void *
      data () const
      {
        return data_;
      }

      hid_t       type;
      size_type   size;
      void        *data_;
    };

    template <typename T>
    struct hdf5_pod_impl : hdf5_pod__
    {
      typedef size_t size_type;
      typedef typename get_type <typename T::lhs_type>::type type_t;

      hdf5_pod_impl (const T &h)
      {
        type  = hdf5::private_::get_hdf5_type_ (type_t ());
        size  = T::static_size;
        data_ = buffer_.data ();

        hdf5::private_::buffer_filler <T::static_size>::fill (h, buffer_);
      }

      boost::array <type_t, T::static_size> buffer_;
    };

    template <typename lhs_t, typename rhs_t>
    struct hdf5_value_holder
    {
      typedef lhs_t lhs_type;
      typedef rhs_t rhs_type;
      
      hdf5_value_holder (const lhs_t &lhs, const rhs_t &rhs)
      : lhs (lhs)
      , rhs (rhs)
      {
      }

      enum { 
        static_size = get_size <lhs_t>::static_size + get_size <rhs_t>::static_size, 
      };

      const lhs_t &lhs;
      const rhs_t &rhs;
    };

    template <typename T>
    struct hdf5_value_holder_unary
    {
      typedef T type;
      hdf5_value_holder_unary (const T &v)
      : v (v)
      {
      }

      enum { static_size = get_size <T>::static_size, };

      const T &v;
    };

    template <typename lhs_t, typename rhs_t>
    struct hdf5_struct_holder
    {
      typedef lhs_t lhs_type;
      typedef rhs_t rhs_type;

      hdf5_struct_holder (const lhs_t &lhs, const rhs_t &rhs)
      : lhs (lhs)
      , rhs (rhs)
      {
      }

      enum { 
        static_size = get_size <lhs_t>::static_size + get_size <rhs_t>::static_size
      };

      typedef typename lhs_type::type_t type_t;

      const lhs_t &lhs;
      const rhs_t &rhs;
    };

  } // namespace private
  } // namespace hdf5

  template <typename T>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const shared_vector <T> &data, 
      hdf5::private_::hdf5_buffer__::size_type stride = 0, 
      hdf5::private_::hdf5_buffer__::size_type offset = 0)
  {
    return hdf5::private_::hdf5_buffer__ (T (), data.data (), data.size (), stride, offset);
  }

  template <typename T>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const seq_vector <T> &data, 
      hdf5::private_::hdf5_buffer__::size_type stride = 0, 
      hdf5::private_::hdf5_buffer__::size_type offset = 0)
  {
    return hdf5::private_::hdf5_buffer__ (T (), data.data (), data.size (), stride, offset);
  }

  template <>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const shared_vector <double> &data,
      hdf5::private_::hdf5_buffer__::size_type stride,
      hdf5::private_::hdf5_buffer__::size_type offset)
  {
    return hdf5::private_::hdf5_buffer__ (float (), data.data (), data.size (), stride, offset);
  }

  template <>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const seq_vector <double> &data,
      hdf5::private_::hdf5_buffer__::size_type stride,
      hdf5::private_::hdf5_buffer__::size_type offset)
  {
    return hdf5::private_::hdf5_buffer__ (float (), data.data (), data.size (), stride, offset);
  }

  template <typename T>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const T &t)
  {
    return hdf5::private_::hdf5_buffer__ (T (), &t, 1, 0, 0);
  }

  template <typename T>
  inline hdf5::private_::hdf5_value_holder_unary <T>
  hdf5_pod (const T &v)
  {
    return hdf5::private_::hdf5_value_holder_unary <T> (v);
  }

  template <typename T>
  struct hdf5_struct
  {
    typedef T type_t;
    enum { static_size = 0, };

    hdf5_struct ()
    : type (get_hdf5_type (T ()))
    , size (0)
    , data_ (0)
    {
    }

    ~hdf5_struct ()
    {
      delete [] static_cast <T *> (data_);
    }

    void
    write (const hdf5::private_::hdf5_value_holder_unary <T> &h)
    {
      using namespace hdf5::private_;
      write_vec (&h.v, 1);
    }

    template <typename L, typename R>
    void
    write (const hdf5::private_::hdf5_value_holder <L, R> &h)
    {
      using namespace hdf5::private_;
      write_pod (hdf5_pod_impl <hdf5_value_holder <L, R> > (h));
    }

    void 
    write (const seq_vector <T> &v)
    {
      write_vec (&v[0], v.size ());
    }

  private:
    struct data_chunk 
    {
      T       *data;
      size_t  count;

      data_chunk (T *data, size_t count)
      : data (data)
      , count (count)
      {
      }
    };

  public:
    void *
    data () const
    {
      if (!data_)
        {
          T *buffer = new T [size];
          for (size_t i = 0, idx = 0, cnt = chunks_.size (); i < cnt; ++i)
            {
              const data_chunk &chunk = chunks_[i];
              for (size_t j = 0; j < chunk.count; ++j, ++idx)
                {
                  buffer[idx] = chunk.data[j];
                }
            }

          data_ = buffer;
        }

      return data_;
    }

  private:

    void
    write_vec (T *t, size_t s)
    {
      T *buffer = new T [s];
      for (size_t i = 0; i < s; ++i)
        buffer[i] = t[i];

      chunks_.push_back (data_chunk (buffer, s));
      size += s;
    }

    template <typename Y>
    void
    write_pod (const Y &pod)
    {
      using namespace hdf5::private_;

      T *buffer = new T [pod.size];
      for (size_t i = 0; i < pod.buffer_.size (); ++i)
        buffer[i] = pod.buffer_[i];

      chunks_.push_back (data_chunk (buffer, pod.size));
      size += pod.size;
    }

  public:
    hid_t                     type;
    size_t                    size;

  private:
    mutable void              *data_;
    seq_vector <data_chunk>   chunks_;
  };

  namespace hdf5 {
  namespace private_ {

    template <size_t s>
    struct hdf5_struct_helper 
    {
      template <typename T>
      static hid_t
      get_type (const T &h)
      {
        return hdf5_struct_helper <T::lhs_type::static_size>::get_type (h.lhs);
      }
      template <typename T>
      static hid_t
      get_size (const T &h)
      {
        return hdf5_struct_helper <T::lhs_type::static_size>::get_size (h.lhs);
      }
      template <typename T, typename B>
      static size_t
      fill_buffer (const T &h, B &buffer, size_t offset)
      {
        offset += hdf5_struct_helper <T::lhs_type::static_size>::fill_buffer (h.lhs, buffer, offset);
        buffer [offset] = h.rhs;
        return offset + 1;
      }
    };

    template <>
    struct hdf5_struct_helper <0>
    {
      template <typename T>
      static hid_t
      get_type (const T &h)
      {
        return h.type;
      }
      template <typename T>
      static hid_t
      get_size (const T &h)
      {
        return h.size;
      }
      template <typename T, typename B>
      static size_t
      fill_buffer (const T &h, B &buffer, size_t offset)
      {
        BS_ASSERT (offset == 0) (offset);
        for (; offset < h.size; ++offset)
          buffer[offset] = static_cast <const typename T::type_t *> (h.data ())[offset];

        return offset;
      }
    };

    struct hdf5_struct__
    {
      template <typename L, typename R>
      hdf5_struct__ (const hdf5_struct_holder <L, R> &h)
      : type (-1)
      , size (get_size <L>::static_size + get_size <R>::static_size)
      , data_ (0)
      {
        typedef hdf5_struct_holder <L, R> T;

        size += hdf5_struct_helper <T::static_size>::get_size (h);
        type  = get_hdf5_type (typename T::type_t ());

        typename T::type_t *buffer = new typename T::type_t [size];
        hdf5_struct_helper <T::static_size>::fill_buffer (h, buffer, 0);

        data_ = buffer;
      }

      void *
      data () const
      {
        return data_;
      }

      hid_t     type;
      size_t    size;
      void      *data_;
    };

  } // namespace private_
  } // namespace hdf5


  struct BS_API_PLUGIN hdf5_file;
  struct BS_API_PLUGIN hdf5_group_v2;
  struct BS_API_PLUGIN hdf5_storage_v2
  {
  private:

    static hdf5_storage_v2 *
    instance ();

    hdf5_storage_v2 ();
    ~hdf5_storage_v2 ();

    struct impl;
    impl *impl_;

    friend struct hdf5_file;
    friend struct hdf5_group_v2;
    friend struct impl;
  };

  struct BS_API_PLUGIN hdf5_file
  {
    hdf5_group_v2 
    operator [] (const std::string &name);

    hdf5_file (const std::string &file_name)
    : file_id_ (-1)
    , file_name_ (file_name)
    {
    }

    ~hdf5_file ();

  private:

    hid_t       file_id_;
    std::string file_name_;

    friend struct hdf5_group_v2;
    friend struct hdf5_property_v2;
    friend struct hdf5_storage_v2;
  };

  struct BS_API_PLUGIN hdf5_group_v2
  {
  private:
    typedef hdf5::private_::hdf5_buffer__ hdf5_buffer_t;
    typedef hdf5::private_::hdf5_pod__    hdf5_pod_t;
    typedef hdf5::private_::hdf5_struct__ hdf5_struct_t;

  public:

    hdf5_group_v2 (const hdf5_file &file, const std::string &name)
    : file_ (file)
    , name_ (name)
    {
    }

    template <typename L>
    hdf5_group_v2
    write (const char *dataset, const hdf5::private_::hdf5_value_holder_unary <L> &h)
    {
      return write_buffer (dataset, hdf5_buffer (h.v));
    }

    template <typename L, typename R>
    hdf5_group_v2 &
    write (const char *dataset, const hdf5::private_::hdf5_value_holder <L, R> &h)
    {
      using namespace hdf5::private_;
      return write_pod (dataset, hdf5_pod_impl <hdf5_value_holder <L, R> > (h));
    }

    hdf5_group_v2 &
    write (const char *dataset, const hdf5_buffer_t &buffer)
    {
      return write_buffer (dataset, buffer);
    }

    template <typename T>
    hdf5_group_v2 &
    write (const char *dataset, const seq_vector <T> &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }

    template <typename T>
    hdf5_group_v2 &
    write (const char *dataset, const hdf5_struct <T> &s)
    {
      return write_buffer (dataset, hdf5_buffer_t (s.data (), s.size, s.type));
    }

    template <typename L, typename R>
    hdf5_group_v2 &
    write (const char *dataset, const hdf5::private_::hdf5_struct_holder <L, R> &s)
    {
      return write_struct (dataset, hdf5_struct_t (s));
    }

  private:

    hdf5_group_v2 &
    write_buffer (const char *dataset, const hdf5_buffer_t &buffer);

    hdf5_group_v2 &
    write_pod (const char *dataset, const hdf5_pod_t &pod);

    hdf5_group_v2 &
    write_struct (const char *dataset, const hdf5_struct_t &s);

  private:
    hdf5_file   file_;
    std::string name_;
  };

  template <typename L, typename R>
  inline hdf5::private_::hdf5_value_holder <L, R>
  operator << (const hdf5::private_::hdf5_value_holder_unary <L> &lhs, const R &rhs)
  {
    return hdf5::private_::hdf5_value_holder <L, R> (lhs.v, rhs);
  }

  template <typename LL, typename LR, typename R>
  inline hdf5::private_::hdf5_value_holder <hdf5::private_::hdf5_value_holder <LL, LR>, R>
  operator << (const hdf5::private_::hdf5_value_holder <LL, LR> &lhs, const R &rhs)
  {
    return hdf5::private_::hdf5_value_holder <hdf5::private_::hdf5_value_holder <LL, LR>, R> (lhs, rhs);
  }

  template <typename T, typename R>
  inline hdf5::private_::hdf5_struct_holder <hdf5_struct <T>, R>
  operator << (const hdf5_struct <T> &lhs, const R &rhs)
  {
    return hdf5::private_::hdf5_struct_holder <hdf5_struct <T>, R> (lhs, rhs);
  }
  template <typename LL, typename LR, typename R>
  inline hdf5::private_::hdf5_struct_holder <hdf5::private_::hdf5_struct_holder <LL, LR>, R>
  operator << (const hdf5::private_::hdf5_struct_holder <LL, LR> &lhs, const R &rhs)
  {
    using namespace hdf5::private_;
    return hdf5_struct_holder <hdf5_struct_holder <LL, LR>, R> (lhs, rhs);
  }

} // namespace blue_sky

#endif
