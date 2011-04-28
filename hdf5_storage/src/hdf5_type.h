/**
 * \file hdf5_type.h
 * \brief Helper to deduce H5:PredType from C++ type
 * \author Vadim Alimguzhin
 * \date 09.10.2008
 */

#ifndef HDF5_TYPE_H_
#define HDF5_TYPE_H_

#include "hdf5.h"
#include <vector>

namespace blue_sky
  {

    enum hdf5_type
    {
      hdf5_char,
      hdf5_int,
      hdf5_uint,
      hdf5_long,
      hdf5_ulong,
      hdf5_float,
      hdf5_double,
    };

  template <typename T>
  struct hdf5_type_helper
  {
    enum { 
      type = hdf5_type_helper <typename T::value_type>::type, 
    };
  };

  template <>
  struct hdf5_type_helper <char>
  {
    enum { type = hdf5_char, };
  };
  template <>
  struct hdf5_type_helper <int>
  {
    enum { type = hdf5_int, };
  };
  template <>
  struct hdf5_type_helper <unsigned int>
  {
    enum { type = hdf5_uint, };
  };
  template <>
  struct hdf5_type_helper <long>
  {
    enum { type = hdf5_long, };
  };
  template <>
  struct hdf5_type_helper <unsigned long>
  {
    enum { type = hdf5_ulong, };
  };
  template <>
  struct hdf5_type_helper <float>
  {
    enum { type = hdf5_float, };
  };
  template <>
  struct hdf5_type_helper <double>
  {
    enum { type = hdf5_double, };
  };

  template<typename T>
  const hid_t & get_hdf5_type()
  {
    BOOST_STATIC_ASSERT (sizeof (T) == sizeof (char));
    return -1;
  }

  template<> inline const hid_t & get_hdf5_type <char> ()
  {
    return H5T_NATIVE_CHAR;
  }
  template<> inline const hid_t & get_hdf5_type <int> ()
  {
    return H5T_NATIVE_INT;
  }
  template<> inline const hid_t & get_hdf5_type <float> ()
  {
    return H5T_NATIVE_FLOAT;
  }
  template<> inline const hid_t & get_hdf5_type <double> ()
  {
    return H5T_NATIVE_DOUBLE;
  }
  template<> inline const hid_t & get_hdf5_type <unsigned long> ()
  {
    return H5T_NATIVE_ULONG;
  }
  template<> inline const hid_t & get_hdf5_type <unsigned int> ()
  {
    return H5T_NATIVE_UINT;
  }
  template <> inline const hid_t & get_hdf5_type <long> ()
  {
    return H5T_NATIVE_LONG;
  }

  template<typename T> inline const hid_t & get_hdf5_type (const T &)
  {
    return get_hdf5_type<T>();
  }
  template<typename T> inline const hid_t & get_hdf5_type (const T x[])
  {
    return get_hdf5_type<T>();
  }
  template<typename T> inline const hid_t & get_hdf5_type (const std::vector <T> &)
  {
    return get_hdf5_type<T>();
  }

} // namespace blue_sky
#endif // #ifndef HDF5_TYPE_H_
