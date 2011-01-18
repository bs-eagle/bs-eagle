/**
 * \file shared_array.h
 * \brief shared array
 * \author Sergey Miryanov
 * \date 28.05.2009
 * */
#ifndef BS_BOS_CORE_DATA_STORAGE_pool_shared_array_H_
#define BS_BOS_CORE_DATA_STORAGE_pool_shared_array_H_

#include <boost/shared_ptr.hpp>
#include "array_ext.h"

namespace blue_sky {

using namespace private_;

  template <typename T>
  struct pool_shared_array 
  {
    // type definitions
    typedef T              value_type;
    typedef T*             iterator;
    typedef const T*       const_iterator;
    typedef T&             reference;
    typedef const T&       const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;

    pool_shared_array (T *e = 0, size_type N = 0)
    : array_ (new array_ext <T> (e, N))
    {
    }

    // iterator support
    iterator begin()
    {
      return array_->begin ();
    }
    const_iterator begin() const
    {
      return array_->begin ();
    }
    iterator end()
    {
      return array_->end ();
    }
    const_iterator end() const
    {
      return array_->end ();
    }

    //reverse_iterator rbegin()
    //{
    //  return reverse_iterator(end());
    //}
    //const_reverse_iterator rbegin() const
    //{
    //  return const_reverse_iterator(end());
    //}
    //reverse_iterator rend()
    //{
    //  return reverse_iterator(begin());
    //}
    //const_reverse_iterator rend() const
    //{
    //  return const_reverse_iterator(begin());
    //}

    // operator[]
    reference operator[](size_type i)
    {
      return array_->operator[] (i);
    }

    const_reference operator[](size_type i) const
    {
      return array_->operator[] (i);
    }

    // at() with range check
    reference at(size_type i)
    {
      return array_->at (i);
    }
    const_reference at(size_type i) const
    {
      return array_->at (i);
    }
    // size is constant
    size_type size() const
    {
      return array_->size ();
    }
    bool empty() const
    {
      return array_->empty ();
    }
    size_type max_size() const
    {
      return array_->max_size ();
    }

    // direct access to data (read-only)
    const T* data() const
    {
      return array_->data ();
    }
    T* data()
    {
      return array_->data ();
    }

    // use array as C array (direct read/write access to data)
    T* c_array()
    {
      return array_->c_array ();
    }
    pool_shared_array <T> &operator= (const pool_shared_array <T> &rhs)
    {
      array_ = rhs.array_;
      return *this;
    }
    // assign one value to all elements
    void assign (const T& value)
    {
      array_->assign (value);
    }

  private:
    boost::shared_ptr <array_ext <T> > array_;
  };


  typedef unsigned char             uchar_t;
  typedef float                     float_t;

  typedef pool_shared_array <uchar_t>    array_uchar_t;
  typedef pool_shared_array <float_t>  array_float_t;

} // namespace blue_sky


#endif // #ifndef BS_BOS_CORE_DATA_STORAGE_pool_shared_array_H_

