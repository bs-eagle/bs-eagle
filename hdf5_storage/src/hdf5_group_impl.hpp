#ifndef HDF5_GROUP_IMPL_HPP_e9a85138_85db_11e0_8944_0ff472cf95ce
#define HDF5_GROUP_IMPL_HPP_e9a85138_85db_11e0_8944_0ff472cf95ce
/**
 *       \file  hdf5_group_impl.hpp
 *      \brief  Impl of hdf5_group_impl_iface
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  24.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "bs_hdf5_storage_v2.h"

namespace blue_sky
{
  struct hdf5_group_impl : hdf5::private_::hdf5_group_impl_iface
  {
    BLUE_SKY_TYPE_DECL (hdf5_group_impl);

    typedef hdf5::private_::hdf5_buffer__ hdf5_buffer_t;
    typedef hdf5::private_::hdf5_pod__    hdf5_pod_t;
    typedef hdf5::private_::hdf5_struct__ hdf5_struct_t;

    virtual void
    init (hdf5_file const &file, std::string const &name);

    virtual void
    write_buffer (const char *dataset, const hdf5_buffer_t &buffer);

    virtual void
    write_pod (const char *dataset, const hdf5_pod_t &pod);

    virtual void
    write_struct (const char *dataset, const hdf5_struct_t &s);

    virtual void
    write_string (const char *dataset, const std::string &value);

    virtual void
    read_buffer (const char *dataset, hdf5_buffer_t &buffer);


    hdf5_file   file_;
    std::string name_;
  };

}

#endif
