/**
 *       \file  hdf5_hid_holder.hpp
 *      \brief  keywords for SCAL
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  28.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "hdf5.h"

namespace blue_sky {
namespace detail {

  enum hid_type {
    hdf5_property_type
    , hdf5_dataspace_type
    , hdf5_dataset_type
    , hdf5_group_type
    , hdf5_datatype_type
  };

  template <hid_type type>
  struct hid_closer
  {
  };

  template <>
  struct hid_closer <hdf5_property_type>
  {
    static void
    close (hid_t v)
    {
      H5Pclose (v);
    }
  };

  template <>
  struct hid_closer <hdf5_dataspace_type>
  {
    static void
    close (hid_t v)
    {
      H5Sclose (v);
    }
  };

  template <>
  struct hid_closer <hdf5_dataset_type>
  {
    static void
    close (hid_t v)
    {
      H5Dclose (v);
    }
  };

  template <>
  struct hid_closer <hdf5_group_type>
  {
    static void
    close (hid_t v)
    {
      H5Gclose (v);
    }
  };

  template <>
  struct hid_closer <hdf5_datatype_type>
  {
    static void
    close (hid_t v)
    {
      H5Tclose (v);
    }
  };

  template <hid_type hid_type_t>
  struct hid_holder
  {
    hid_holder (hid_t h)
    : h_ (h)
    {
    }

    ~hid_holder ()
    {
      if (h_ >= 0)
        hid_closer <hid_type_t>::close (h_);
    }

    hid_holder &
    operator= (hid_t h)
    {
      if (h_ >= 0)
        hid_closer <hid_type_t>::close (h_);

      h_ = h;
      return *this;
    }

    operator hid_t () 
    {
      return h_;
    }

    bool
    valid () const
    {
      return h_ >= 0;
    }

    hid_t h_;
  };

} // namespace detail

typedef detail::hid_holder <detail::hdf5_group_type>      hid_group_t;
typedef detail::hid_holder <detail::hdf5_dataspace_type>  hid_dspace_t;
typedef detail::hid_holder <detail::hdf5_dataset_type>    hid_dset_t;
typedef detail::hid_holder <detail::hdf5_property_type>   hid_property_t;
typedef detail::hid_holder <detail::hdf5_datatype_type>   hid_dtype_t;

} // namespace blue_sky

