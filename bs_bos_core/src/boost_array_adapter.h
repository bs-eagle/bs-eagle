/**
 *       \file  boost_array_adapter.h
 *      \brief  Converts calc_model_data members (boost::array) to
 *              shared_vector. This file should be included first.
 *     \author  Sergey Miryanov
 *       \date  02.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_EAGLE_BOOST_ARRAY_ADAPTER_H_
#define BS_EAGLE_BOOST_ARRAY_ADAPTER_H_

#include <boost/array.hpp>
#include <boost/mpl/vector.hpp>

#include "shared_vector.h"
#include "constants.h"
#include "strategies.h"
#include "calc_model_data.h"

namespace blue_sky {
namespace detail {

  template <typename T, typename Y, size_t N>
  struct boost_array_adapter__
  {
    boost_array_adapter__ (boost::array <Y, N> T::*array)
    : array (array)
    {
    }

    shared_vector <Y>
    operator () (T *t)
    {
      return shared_array (((*t).*array).data (), ((*t).*array).size ());
    }

    boost::array <Y, N> T::*array;
  };

  template <typename T, typename Y, size_t N>
  boost_array_adapter__ <T, Y, N> 
  boost_array_adapter (boost::array <Y, N> T::*array)
  {
    return boost_array_adapter__ <T, Y, N> (array);
  }

} // namespace detail
} // namespace blue_sky

#ifdef BSPY_EXPORTING_PLUGIN
namespace boost {
namespace python {
namespace detail {

  template <typename S, typename T, size_t N>
  mpl::vector <blue_sky::shared_vector <T>, blue_sky::calc_model_data <S> *>
  get_signature (const blue_sky::detail::boost_array_adapter__ <blue_sky::calc_model_data <S>, T, N> &, 
    blue_sky::calc_model_data <S> *)
  {
    return mpl::vector <blue_sky::shared_vector <T>, blue_sky::calc_model_data <S> *> ();
  }

} // namespace detail
} // namespace python
} // namespace boost

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#endif

#endif // #ifndef BS_EAGLE_BOOST_ARRAY_ADAPTER_H_

