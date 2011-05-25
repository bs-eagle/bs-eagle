#ifndef HDF5_STORAGE_V2_IMPL_HPP_32ed5d72_85e4_11e0_9665_6fa901f622fd
#define HDF5_STORAGE_V2_IMPL_HPP_32ed5d72_85e4_11e0_9665_6fa901f622fd
/**
 *       \file  hdf5_storage_v2_impl.hpp
 *      \brief  Some bs_hdf5_storage_v2 implementation defaults
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  24.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include <boost/type_traits.hpp>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include "auto_value.h"

#include "hdf5_type.h"

namespace blue_sky {

  struct BS_API_PLUGIN hdf5_file;
  struct BS_API_PLUGIN hdf5_group_v2;

  namespace hdf5 {
  namespace private_ {

    struct hdf5_buffer__ 
    {
      typedef size_t size_type;

      template <typename T, typename Y>
      hdf5_buffer__ (const Y &, const T *data, size_type size_, size_type stride, size_type offset)
      : type (hdf5_type_helper <Y>::type)
      , size (stride ? (size_ / stride) : size_)
      , owner_ (stride || !boost::is_same <typename boost::remove_cv <T>::type, typename boost::remove_cv <Y>::type> ())
      , data_ (!owner_ ? (void *)data : (void *)(new Y[ size]))
      {
        if (owner_)
          {
            if (!stride)
              stride = 1;

            for (size_type i = offset, j = 0; i < size_; i += stride, ++j)
              {
                static_cast <Y *> (data_)[j] = static_cast <Y> (data[i]);
              }
          }
      }

      hdf5_buffer__ (void *data, size_type size, hdf5_type type)
      : type (type)
      , size (size)
      , owner_ (false)
      , data_ (data)
      {
      }

      ~hdf5_buffer__ ()
      {
        if (owner_)
          delete [] reinterpret_cast <char *> (data_);
      }

      void *
      data () const
      {
        return data_;
      }

    public:
      hdf5_type   type;
      size_type   size;

      bool        owner_;
      void        *data_;
    };

    template <typename T>
    struct get_size
    {
      enum { static_size = T::static_size, };
    };

    template <>
    struct get_size <size_t>
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <int>
    {
      enum { static_size = 1,  };
    };

    template <>
    struct get_size <float>
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <double>
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <auto_value <double> >
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <auto_value <float> >
    {
      enum { static_size = 1, };
    };
    template <>
    struct get_size <long>
    {
      enum { static_size = 1, };
    };

    template <typename T>
    struct get_type
    {
      typedef typename get_type <typename T::lhs_type>::type type;
    };

    template <>
    struct get_type <size_t>
    {
      typedef size_t type;
    };
    template <>
    struct get_type <int>
    {
      typedef int type;
    };
    template <>
    struct get_type <float>
    {
      typedef float type;
    };
    template <>
    struct get_type <double>
    {
      typedef double type;
    };
    template <>
    struct get_type <auto_value <double> >
    {
      typedef double type;
    };
    template <>
    struct get_type <auto_value <float> >
    {
      typedef float type;
    };
    template <>
    struct get_type <long>
    {
      typedef long type;
    };

    template <size_t s>
    struct buffer_filler 
    {
      template <typename R, typename B>
      static void
      fill (const R &h, B &buffer)
      {
        buffer_filler <R::lhs_type::static_size>::fill (h.lhs, buffer);
        buffer[R::static_size - 1] = h.rhs;
      }
    };
    template <>
    struct buffer_filler <2>
    {
      template <typename R, typename B>
      static void
      fill (const R &h, B &buffer)
      {
        buffer[0] = h.lhs;
        buffer[1] = h.rhs;
      }
    };

    struct hdf5_pod__
    {
      typedef size_t size_type;

      hdf5_pod__ ()
      : type (hdf5_none)
      , size (0)
      , data_ (0)
      {
      }

      void *
      data () const
      {
        return data_;
      }

      hdf5_type   type;
      size_type   size;
      void        *data_;
    };

    template <typename T>
    struct hdf5_pod_impl : hdf5_pod__
    {
      typedef size_t size_type;
      typedef typename get_type <typename T::lhs_type>::type type_t;

      hdf5_pod_impl (const T &h)
      {
        type  = hdf5_type_helper <type_t>::type;
        size  = T::static_size;
        data_ = buffer_.data ();

        hdf5::private_::buffer_filler <T::static_size>::fill (h, buffer_);
      }

      boost::array <type_t, T::static_size> buffer_;
    };

    template <typename lhs_t, typename rhs_t>
    struct hdf5_value_holder
    {
      typedef lhs_t lhs_type;
      typedef rhs_t rhs_type;
      
      hdf5_value_holder (const lhs_t &lhs, const rhs_t &rhs)
      : lhs (lhs)
      , rhs (rhs)
      {
      }

      enum { 
        static_size = get_size <lhs_t>::static_size + get_size <rhs_t>::static_size, 
      };

      const lhs_t &lhs;
      const rhs_t &rhs;
    };

    template <typename T>
    struct hdf5_value_holder_unary
    {
      typedef T type;
      hdf5_value_holder_unary (const T &v)
      : v (v)
      {
      }

      enum { static_size = get_size <T>::static_size, };

      const T &v;
    };

    template <typename lhs_t, typename rhs_t>
    struct hdf5_struct_holder
    {
      typedef lhs_t lhs_type;
      typedef rhs_t rhs_type;

      hdf5_struct_holder (const lhs_t &lhs, const rhs_t &rhs)
      : lhs (lhs)
      , rhs (rhs)
      {
      }

      enum { 
        static_size = get_size <lhs_t>::static_size + get_size <rhs_t>::static_size
      };

      typedef typename lhs_type::type_t type_t;

      const lhs_t &lhs;
      const rhs_t &rhs;
    };

  } // namespace private
  } // namespace hdf5

  template <typename T>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const shared_vector <T> &data, 
      hdf5::private_::hdf5_buffer__::size_type stride = 0, 
      hdf5::private_::hdf5_buffer__::size_type offset = 0)
  {
    return hdf5::private_::hdf5_buffer__ (T (), data.data (), data.size (), stride, offset);
  }

  template <>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const shared_vector <double> &data,
      hdf5::private_::hdf5_buffer__::size_type stride,
      hdf5::private_::hdf5_buffer__::size_type offset)
  {
    return hdf5::private_::hdf5_buffer__ (float (), data.data (), data.size (), stride, offset);
  }
  template <>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const shared_vector <main_var_type> &data,
    hdf5::private_::hdf5_buffer__::size_type stride,
    hdf5::private_::hdf5_buffer__::size_type offset)
  {
    return hdf5::private_::hdf5_buffer__ (int (), data.data (), data.size (), stride, offset);
  }

  template <typename T>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const T &t)
  {
    return hdf5::private_::hdf5_buffer__ (T (), &t, 1, 0, 0);
  }

  template <>
  inline hdf5::private_::hdf5_buffer__ 
  hdf5_buffer (const spv_double &data)
  {
    return hdf5::private_::hdf5_buffer__ (float (), data->data (), data->size (), 0, 0);
  }
#if !T_FLOAT_IS_DOUBLE 
  template <>
  inline hdf5::private_::hdf5_buffer__ 
  hdf5_buffer (const spv_float &data)
  {
    return hdf5::private_::hdf5_buffer__ (float (), data->data (), data->size (), 0, 0);
  }
#endif
  template <>
  inline hdf5::private_::hdf5_buffer__ 
  hdf5_buffer (const spv_long &data)
  {
    return hdf5::private_::hdf5_buffer__ (long (), data->data (), data->size (), 0, 0);
  }
  template <>
  inline hdf5::private_::hdf5_buffer__ 
  hdf5_buffer (const spv_int &data)
  {
    return hdf5::private_::hdf5_buffer__ (int (), data->data (), data->size (), 0, 0);
  }
  template <typename T>
  inline hdf5::private_::hdf5_buffer__ 
  hdf5_buffer (const std::vector <T> &data)
  {
    return hdf5::private_::hdf5_buffer__ (T (), &(data[0]), data.size (), 0, 0);
  }
  template <>
  inline hdf5::private_::hdf5_buffer__ 
  hdf5_buffer (const std::vector <main_var_type> &data)
  {
    return hdf5::private_::hdf5_buffer__ (int (), &(data[0]), data.size (), 0, 0);
  }

  template <typename T, size_t N>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (boost::array <T, N> const &data)
  {
    return hdf5::private_::hdf5_buffer__ (T (), data.data (), N, 0, 0);
  }


  template <typename T>
  inline hdf5::private_::hdf5_value_holder_unary <T>
  hdf5_pod (const T &v)
  {
    return hdf5::private_::hdf5_value_holder_unary <T> (v);
  }

  template <typename T>
  struct hdf5_struct
  {
    typedef T type_t;
    enum { static_size = 0, };

    hdf5_struct ()
    : type (hdf5_type_helper <T>::type)
    , size (0)
    , data_ (0)
    {
    }

    ~hdf5_struct ()
    {
      delete [] static_cast <T *> (data_);
    }

    void
    write (const hdf5::private_::hdf5_value_holder_unary <T> &h)
    {
      using namespace hdf5::private_;
      write_vec (&h.v, 1);
    }

    template <typename L, typename R>
    void
    write (const hdf5::private_::hdf5_value_holder <L, R> &h)
    {
      using namespace hdf5::private_;
      write_pod (hdf5_pod_impl <hdf5_value_holder <L, R> > (h));
    }

    void 
    write (const shared_vector <T> &v)
    {
      write_vec (&v[0], v.size ());
    }

  private:
    struct data_chunk 
    {
      T       *data;
      size_t  count;

      data_chunk (T *data, size_t count)
      : data (data)
      , count (count)
      {
      }
    };

  public:
    void *
    data () const
    {
      if (!data_)
        {
          T *buffer = new T [size];
          for (size_t i = 0, idx = 0, cnt = chunks_.size (); i < cnt; ++i)
            {
              const data_chunk &chunk = chunks_[i];
              for (size_t j = 0; j < chunk.count; ++j, ++idx)
                {
                  buffer[idx] = chunk.data[j];
                }
            }

          data_ = buffer;
        }

      return data_;
    }

  private:

    void
    write_vec (T *t, size_t s)
    {
      T *buffer = new T [s];
      for (size_t i = 0; i < s; ++i)
        buffer[i] = t[i];

      chunks_.push_back (data_chunk (buffer, s));
      size += s;
    }

    template <typename Y>
    void
    write_pod (const Y &pod)
    {
      using namespace hdf5::private_;

      T *buffer = new T [pod.size];
      for (size_t i = 0; i < pod.buffer_.size (); ++i)
        buffer[i] = pod.buffer_[i];

      chunks_.push_back (data_chunk (buffer, pod.size));
      size += pod.size;
    }

  public:
    hdf5_type                   type;
    size_t                      size;

  private:
    mutable void                *data_;
    shared_vector <data_chunk>  chunks_;
  };

  namespace hdf5 {
  namespace private_ {

    template <size_t s>
    struct hdf5_struct_helper 
    {
      template <typename T>
      static hdf5_type
      get_type (const T &h)
      {
        return hdf5_struct_helper <T::lhs_type::static_size>::get_type (h.lhs);
      }
      template <typename T>
      static size_t
      get_size (const T &h)
      {
        return hdf5_struct_helper <T::lhs_type::static_size>::get_size (h.lhs);
      }
      template <typename T, typename B>
      static size_t
      fill_buffer (const T &h, B &buffer, size_t offset)
      {
        offset += hdf5_struct_helper <T::lhs_type::static_size>::fill_buffer (h.lhs, buffer, offset);
        buffer [offset] = h.rhs;
        return offset + 1;
      }
    };

    template <>
    struct hdf5_struct_helper <0>
    {
      template <typename T>
      static hdf5_type
      get_type (const T &h)
      {
        return h.type;
      }
      template <typename T>
      static size_t
      get_size (const T &h)
      {
        return h.size;
      }
      template <typename T, typename B>
      static size_t
      fill_buffer (const T &h, B &buffer, size_t offset)
      {
        BS_ASSERT (offset == 0) (offset);
        for (; offset < h.size; ++offset)
          buffer[offset] = static_cast <const typename T::type_t *> (h.data ())[offset];

        return offset;
      }
    };

    struct hdf5_struct__
    {
      template <typename L, typename R>
      hdf5_struct__ (const hdf5_struct_holder <L, R> &h)
      : type (hdf5_none)
      , size (get_size <L>::static_size + get_size <R>::static_size)
      , data_ (0)
      {
        typedef hdf5_struct_holder <L, R> T;

        size += hdf5_struct_helper <T::static_size>::get_size (h);
        type  = hdf5_type_helper <typename T::type_t>::type;

        typename T::type_t *buffer = new typename T::type_t [size];
        hdf5_struct_helper <T::static_size>::fill_buffer (h, buffer, 0);

        data_ = buffer;
      }

      void *
      data () const
      {
        return data_;
      }

      hdf5_type type;
      size_t    size;
      void      *data_;
    };

    struct hdf5_group_impl_iface : objbase
    {
      virtual ~hdf5_group_impl_iface () {}

      virtual void
      init (hdf5_file const &file, std::string const &name) = 0;

      virtual void
      write_buffer (const char *dataset, hdf5_buffer__ const &buffer) = 0;

      virtual void
      write_pod (const char *dataset, hdf5_pod__ const &pod) = 0;

      virtual void
      write_struct (const char *dataset, hdf5_struct__ const &s) = 0;

      virtual void
      write_string (const char *dataset, std::string const &v) = 0;
    };

  } // namespace private_
  } // namespace hdf5

}

#endif
