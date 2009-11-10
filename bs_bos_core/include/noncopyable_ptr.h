/**
 * \file noncopyable_ptr.h
 * \brief type that act as a pointer but disable copying of owned ptr
 * \author Sergey Miryanov
 * \date 12.02.2009
 * */
#ifndef BS_NONCOPYABLE_PTR_H_
#define BS_NONCOPYABLE_PTR_H_

namespace blue_sky {

  template <typename T>
  struct nc_ptr : boost::noncopyable
  {

    explicit nc_ptr (T *ptr)
    : ptr_ (ptr)
    {
    }

    T *operator-> () 
    {
      return ptr_;
    }
    const T *operator-> () const
    {
      return ptr_;
    }

    bool
    operator ! () const
    {
      return !ptr_;
    }

    T *
    get ()
    {
      return ptr_;
    }

  private:

    T *ptr_;
  };


} // namespace blue_sky

#endif  // BS_NONCOPYABLE_PTR_H_

