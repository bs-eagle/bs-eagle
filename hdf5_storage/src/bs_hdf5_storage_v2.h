/**
 * */
#ifndef BS_EAGLE_HDF5_STORAGE_V2_H_
#define BS_EAGLE_HDF5_STORAGE_V2_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "bs_kernel.h"
#include "shared_vector.h"
#include "constants.h"
#include "conf.h"

#include "hdf5_storage_v2_impl.hpp"

namespace blue_sky {

  struct BS_API_PLUGIN hdf5_group_v2
  {
  private:
    typedef hdf5::private_::hdf5_buffer__ hdf5_buffer_t;
    typedef hdf5::private_::hdf5_pod__    hdf5_pod_t;
    typedef hdf5::private_::hdf5_struct__ hdf5_struct_t;
    typedef hdf5::private_::hdf5_group_impl_iface impl_t;

  public:

    hdf5_group_v2 (const hdf5_file &file, const std::string &name)
    : impl_ (BS_KERNEL.create_object ("hdf5_group_impl"))
    {
      BS_ASSERT (impl_);
      impl_->init (file, name);
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
    write (const char *dataset, const shared_vector <T> &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }
    hdf5_group_v2 &
    write (const char *dataset, const spv_double &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }
#if !T_FLOAT_IS_DOUBLE
    hdf5_group_v2 &
    write (const char *dataset, const spv_float &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }
#endif
    hdf5_group_v2 &
    write (const char *dataset, const spv_long &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }
    hdf5_group_v2 &
    write (const char *dataset, const spv_int &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }
    template <typename T>
    hdf5_group_v2 &
    write (const char *dataset, const std::vector <T> &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }
    template <typename T, size_t N>
    hdf5_group_v2 &
    write (const char *dataset, boost::array <T, N> const &data)
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

    hdf5_group_v2 &
    write (const char *dataset, const std::string &value)
    {
      return write_string (dataset, value);
    }

  private:

    hdf5_group_v2 &
    write_buffer (const char *dataset, const hdf5_buffer_t &buffer)
    {
      impl_->write_buffer (dataset, buffer);
      return *this;
    }

    hdf5_group_v2 &
    write_pod (const char *dataset, const hdf5_pod_t &pod)
    {
      impl_->write_pod (dataset, pod);
      return *this;
    }

    hdf5_group_v2 &
    write_struct (const char *dataset, const hdf5_struct_t &s)
    {
      impl_->write_struct (dataset, s);
      return *this;
    }

    hdf5_group_v2 &
    write_string (const char *dataset, const std::string &value)
    {
      impl_->write_string (dataset, value);
      return *this;
    }

  private:
    smart_ptr <impl_t> impl_;
  };


  struct BS_API_PLUGIN hdf5_file
  {
    hdf5_group_v2 
    operator [] (const std::string &name)
    {
      return hdf5_group_v2 (*this, name);
    }

    hdf5_file (const std::string &file_name)
    : file_id_ (-1)
    , file_name_ (file_name)
    {
    }

  private:

    int         file_id_;
    std::string file_name_;

    friend struct hdf5_group_v2;
    friend struct hdf5_storage_v2;
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

#endif // #ifndef BS_EAGLE_HDF5_STORAGE_V2_H_

