#ifndef HDF5_TYPE_TO_HID_HPP_9d584bb4_85ca_11e0_889a_4b717eb768b0
#define HDF5_TYPE_TO_HID_HPP_9d584bb4_85ca_11e0_889a_4b717eb768b0
/**
 *       \file  hdf5_type_to_hid.hpp
 *      \brief  Helper to convert hdf5_type to hid_t
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  24.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "hdf5_type.h"
#include "hdf5.h"

#include "throw_exception.h"

namespace blue_sky {

  template <hdf5_type T>
  const hid_t & get_hdf5_type()
  {
    BOOST_STATIC_ASSERT (T == hdf5_none);
    return -1;
  }

  template<> inline const hid_t & get_hdf5_type <hdf5_char> ()
  {
    return H5T_NATIVE_CHAR;
  }
  template<> inline const hid_t & get_hdf5_type <hdf5_int> ()
  {
    return H5T_NATIVE_INT;
  }
  template<> inline const hid_t & get_hdf5_type <hdf5_float> ()
  {
    return H5T_NATIVE_FLOAT;
  }
  template<> inline const hid_t & get_hdf5_type <hdf5_double> ()
  {
    return H5T_NATIVE_DOUBLE;
  }
  template<> inline const hid_t & get_hdf5_type <hdf5_ulong> ()
  {
    return H5T_NATIVE_ULONG;
  }
  template<> inline const hid_t & get_hdf5_type <hdf5_uint> ()
  {
    return H5T_NATIVE_UINT;
  }
  template <> inline const hid_t & get_hdf5_type <hdf5_long> ()
  {
    return H5T_NATIVE_LONG;
  }

  inline const hid_t & 
  get_hdf5_type (hdf5_type const &type)
  {
    switch (type)
      {
      case hdf5_char  :   return get_hdf5_type <hdf5_char> ();
      case hdf5_int   :   return get_hdf5_type <hdf5_int> ();
      case hdf5_uint  :   return get_hdf5_type <hdf5_uint> ();
      case hdf5_long  :   return get_hdf5_type <hdf5_long> ();
      case hdf5_ulong :   return get_hdf5_type <hdf5_ulong> ();
      case hdf5_float :   return get_hdf5_type <hdf5_float> ();
      case hdf5_double:   return get_hdf5_type <hdf5_double> ();
      default :
        bs_throw_exception (boost::format ("unknown hdf5_type: %d") % type);
      };
  }
}

#endif
