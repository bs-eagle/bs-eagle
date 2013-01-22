/**
 * \file hdf5_type.h
 * \brief Helper to deduce H5:PredType from C++ type
 * \author Vadim Alimguzhin
 * \date 09.10.2008
 */

#ifndef HDF5_TYPE_H_
#define HDF5_TYPE_H_

//#include <vector>

namespace blue_sky
  {

    enum hdf5_type
    {
      hdf5_none,
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
    static const hdf5_type type = hdf5_type_helper <typename T::value_type>::type;
  };

  //template< >
  template< class T >
  struct hdf5_type_helper< smart_ptr< T, true > > {
    static const hdf5_type type = hdf5_type_helper <typename T::value_type>::type;
  };

  template <>
  struct hdf5_type_helper <char>
  {
    static const hdf5_type type = hdf5_char;
  };
  template <>
  struct hdf5_type_helper <int>
  {
    static const hdf5_type type = hdf5_int;
  };
  template <>
  struct hdf5_type_helper <unsigned int>
  {
    static const hdf5_type type = hdf5_uint;
  };
  template <>
  struct hdf5_type_helper <long>
  {
    static const hdf5_type type = hdf5_long;
  };
  template <>
  struct hdf5_type_helper <unsigned long>
  {
    static const hdf5_type type = hdf5_ulong;
  };
  template <>
  struct hdf5_type_helper <float>
  {
    static const hdf5_type type = hdf5_float;
  };
  template <>
  struct hdf5_type_helper <double>
  {
    static const hdf5_type type = hdf5_double;
  };

  //template<typename T> inline const hid_t & get_hdf5_type (const T &)
  //{
  //  return get_hdf5_type<T>();
  //}
  //template<typename T> inline const hid_t & get_hdf5_type (const T x[])
  //{
  //  return get_hdf5_type<T>();
  //}
  //template<typename T> inline const hid_t & get_hdf5_type (const std::vector <T> &)
  //{
  //  return get_hdf5_type<T>();
  //}

} // namespace blue_sky
#endif // #ifndef HDF5_TYPE_H_
