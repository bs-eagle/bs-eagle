/**
 * \file
 * \brief spec for vectors assign
 * \author Sergey Miryanov
 * \date 02.03.2009
 * */
#ifndef BS_VECTOR_ASSIGN_H_
#define BS_VECTOR_ASSIGN_H_

#include "force_inline.h"

namespace blue_sky {

  template <typename item_t, size_t size, typename val_t>
  BS_FORCE_INLINE void
  assign (boost::array <item_t, size> &array, const val_t &val)
  {
    memset (&array[0], val, sizeof (item_t) * size);
  }

  template <typename item_t, typename val_t>
  BS_FORCE_INLINE void
  assign (shared_vector <item_t> &array, const val_t &val)
  {
    size_t size = array.size ();
    if (size)
      memset (&array[0], val, sizeof (item_t) * array.size ());
  }

  template <typename item_t, typename val_t, typename size_type>
  BS_FORCE_INLINE void
  assign (shared_vector <item_t> &array, size_type size, const val_t &val)
  {
    array.resize (size);
    assign (array, val);
  }

} // namespace blue_sky

#endif  // #ifndef BS_VECTOR_ASSIGN_H_

