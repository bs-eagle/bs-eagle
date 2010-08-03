/**
 * \file array_ext.h
 * \brief extended boost::array
 * \author Nikonov Maxim
 * \modification Sergey Miryanov
 * \date 20.02.2009
 * */

#ifndef BOS_ARRAY_EXT_H
#define BOS_ARRAY_EXT_H

namespace blue_sky {

  template<class T>
  class array_ext
  {
  public:
    // type definitions
    typedef T              value_type;
    typedef T*             iterator;
    typedef const T*       const_iterator;
    typedef T&             reference;
    typedef const T&       const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;

    // iterator support
    iterator begin()
    {
      return elems;
    }
    const_iterator begin() const
    {
      return elems;
    }
    iterator end()
    {
      return elems+N;
    }
    const_iterator end() const
    {
      return elems+N;
    }

    // reverse iterator support
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    array_ext (T *e = 0, size_t N = 0) 
    : elems (e)
    , N (N)
    {
    }

    T *
    init (T *d)
    {
      elems = d;
      return elems + N;
    }

    reverse_iterator rbegin()
    {
      return reverse_iterator(end());
    }
    const_reverse_iterator rbegin() const
    {
      return const_reverse_iterator(end());
    }
    reverse_iterator rend()
    {
      return reverse_iterator(begin());
    }
    const_reverse_iterator rend() const
    {
      return const_reverse_iterator(begin());
    }

    // operator[]
    reference operator[](size_type i)
    {
      BS_ASSERT( i < N && N != 0 && "out of range" );
      return elems[i];
    }

    const_reference operator[](size_type i) const
    {
      BS_ASSERT( i < N && N != 0 && "out of range" );
      return elems[i];
    }

    // at() with range check
    reference at(size_type i)
    {
      rangecheck(i);
      return elems[i];
    }
    const_reference at(size_type i) const
    {
      rangecheck(i);
      return elems[i];
    }

    // front() and back()
    reference front()
    {
      return elems[0];
    }

    const_reference front() const
    {
      return elems[0];
    }

    reference back()
    {
      return elems[N-1];
    }

    const_reference back() const
    {
      return elems[N-1];
    }

    // size is constant
    size_type size() const
    {
      return N;
    }
    bool empty() const
    {
      return N == 0;
    }
    size_type max_size() const
    {
      return N;
    }

    // swap (note: linear complexity)
    void swap (array_ext<T> &y)
    {
      BS_ASSERT (y.size () == size ()) (y.size ()) (size ());
      std::swap_ranges(begin(),end(), y.begin());
    }

    // direct access to data (read-only)
    const T* data() const
    {
      return elems;
    }
    T* data()
    {
      return elems;
    }

    // use array as C array (direct read/write access to data)
    T* c_array()
    {
      return elems;
    }

    //// assignment with type conversion
    //template <typename T2>
    //array_ext<T> &operator= (const array_ext<T2>& rhs)
    //{
    //  BS_ASSERT (rhs.size () == size ()) (rhs.size ()) (size ());
    //  std::copy(rhs.begin(),rhs.end(), begin());
    //  return *this;
    //}

    array_ext <T> &operator= (array_ext <T> &rhs)
    {
      BS_ASSERT (rhs.size () == size ()) (rhs.size ()) (size ());
      elems = rhs.elems;
      N = rhs.N;
      return *this;
    }

    // assign one value to all elements
    void assign (const T& value)
    {
      std::fill_n(begin(),size(),value);
    }

    // check range (may be private because it is static)
    void rangecheck (size_type i)
    {
      if (i >= size())
        {
          throw std::out_of_range("array<>: index out of range");
        }
    }
  private:
    T             *elems;
    size_t        N;
  };

  //template< class T >
  //class array_ext< T, 0 >
  //  {

  //  public:
  //    // type definitions
  //    typedef T              value_type;
  //    typedef T*             iterator;
  //    typedef const T*       const_iterator;
  //    typedef T&             reference;
  //    typedef const T&       const_reference;
  //    typedef std::size_t    size_type;
  //    typedef std::ptrdiff_t difference_type;

  //    // iterator support
  //    iterator begin()
  //    {
  //      return iterator( reinterpret_cast< T * >( this ) );
  //    }
  //    const_iterator begin() const
  //      {
  //        return const_iterator(  reinterpret_cast< const T * >( this ) );
  //      }
  //    iterator end()
  //    {
  //      return begin();
  //    }
  //    const_iterator end() const
  //      {
  //        return begin();
  //      }

  //    // reverse iterator support
  //    typedef std::reverse_iterator<iterator> reverse_iterator;
  //    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  //    reverse_iterator rbegin()
  //    {
  //      return reverse_iterator(end());
  //    }
  //    const_reverse_iterator rbegin() const
  //      {
  //        return const_reverse_iterator(end());
  //      }
  //    reverse_iterator rend()
  //    {
  //      return reverse_iterator(begin());
  //    }
  //    const_reverse_iterator rend() const
  //      {
  //        return const_reverse_iterator(begin());
  //      }

  //    // operator[]
  //    reference operator[](size_type i)
  //    {
  //      return failed_rangecheck();
  //    }

  //    const_reference operator[](size_type i) const
  //      {
  //        return failed_rangecheck();
  //      }

  //    // at() with range check
  //    reference at(size_type i)
  //    {
  //      return failed_rangecheck();
  //    }
  //    const_reference at(size_type i) const
  //      {
  //        return failed_rangecheck();
  //      }

  //    // front() and back()
  //    reference front()
  //    {
  //      return failed_rangecheck();
  //    }

  //    const_reference front() const
  //      {
  //        return failed_rangecheck();
  //      }

  //    reference back()
  //    {
  //      return failed_rangecheck();
  //    }

  //    const_reference back() const
  //      {
  //        return failed_rangecheck();
  //      }

  //    // size is constant
  //    static size_type size()
  //    {
  //      return 0;
  //    }
  //    static bool empty()
  //    {
  //      return true;
  //    }
  //    static size_type max_size()
  //    {
  //      return 0;
  //    }
  //    enum { static_size = 0 };

  //    void swap (array_ext<T,0>& y)
  //    {
  //    }

  //    // direct access to data (read-only)
  //    const T* data() const
  //      {
  //        return 0;
  //      }
  //    T* data()
  //    {
  //      return 0;
  //    }

  //    // use array as C array (direct read/write access to data)
  //    T* c_array()
  //    {
  //      return 0;
  //    }

  //    // assignment with type conversion
  //    template <typename T2>
  //    array_ext<T,0>& operator= (const array_ext<T2,0>& )
  //    {
  //      return *this;
  //    }

  //    // assign one value to all elements
  //    void assign (const T& )
  //    {   }

  //    // check range (may be private because it is static)
  //    static reference failed_rangecheck ()
  //    {
  //      std::out_of_range e("attempt to access element of an empty array");
  //      boost::throw_exception(e);
  //      //
  //      // We need to return something here to keep
  //      // some compilers happy: however we will never
  //      // actually get here....
  //      //
  //      static T placeholder;
  //      return placeholder;
  //    }
  //  };

  //// comparisons
  //template<class T, std::size_t N>
  //bool operator== (const array_ext<T,N>& x, const array_ext<T,N>& y)
  //{
  //  return std::equal(x.begin(), x.end(), y.begin());
  //}
  //template<class T, std::size_t N>
  //bool operator< (const array_ext<T,N>& x, const array_ext<T,N>& y)
  //{
  //  return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
  //}
  //template<class T, std::size_t N>
  //bool operator!= (const array_ext<T,N>& x, const array_ext<T,N>& y)
  //{
  //  return !(x==y);
  //}
  //template<class T, std::size_t N>
  //bool operator> (const array_ext<T,N>& x, const array_ext<T,N>& y)
  //{
  //  return y<x;
  //}
  //template<class T, std::size_t N>
  //bool operator<= (const array_ext<T,N>& x, const array_ext<T,N>& y)
  //{
  //  return !(y<x);
  //}
  //template<class T, std::size_t N>
  //bool operator>= (const array_ext<T,N>& x, const array_ext<T,N>& y)
  //{
  //  return !(x<y);
  //}

  // global swap()
  template<class T>
  inline void swap (array_ext<T>& x, array_ext<T>& y)
  {
    BS_ASSERT (x.size () == y.size ()) (x.size ()) (y.size ());
    x.swap(y);
  }
}

#endif // BOS_ARRAY_EXT_H
