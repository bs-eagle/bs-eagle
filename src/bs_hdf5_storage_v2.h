/**
 * */
#ifndef BS_EAGLE_HDF5_STORAGE_V2_H_
#define BS_EAGLE_HDF5_STORAGE_V2_H_

#include "shared_vector.h"
#include "constants.h"

#include <boost/type_traits.hpp>
#include <boost/array.hpp>
#include "auto_value.h"

#include "bs_hdf5_storage.h"

namespace blue_sky {

  namespace hdf5 {
  namespace private_ {

    template <typename T>
    const hid_t &get_hdf5_type_ (const T &)
    {
      return get_hdf5_type <T> ();
    }
    template <>
    inline const hid_t &
    get_hdf5_type_ <main_var_type> (const main_var_type &)
    {
      return H5T_NATIVE_UINT;
    }

    template <typename T>
    inline const hid_t &
    get_hdf5_type_ (const shared_vector <T> &)
    {
      return get_hdf5_type_ <T> (T ());
    }


    struct hdf5_buffer__ 
    {
      typedef size_t size_type;

      template <typename T, typename Y>
      hdf5_buffer__ (const Y &, const T *data, size_type size, size_type stride, size_type offset)
      : type (get_hdf5_type_ (Y ()))
      , size (size)
      , owner_ (stride || !boost::is_same <typename boost::remove_cv <T>::type, typename boost::remove_cv <Y>::type> ())
      , data_ (!owner_ ? (void *)data : (void *)(new Y[stride ? (size / stride) : size]))
      {
        if (owner_)
          {
            if (!stride)
              stride = 1;

            size = size / stride;
            for (size_type i = 0, j = 0; i < size; i += stride, ++j)
              {
                static_cast <Y *> (data_)[j] = static_cast <Y> (data[i]);
              }
          }
      }

      hdf5_buffer__ (void *data, size_type size, hid_t type)
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
      hid_t       type;
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
      : type (0)
      , size (0)
      , data_ (0)
      {
      }

      void *
      data () const
      {
        return data_;
      }

      hid_t       type;
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
        type  = hdf5::private_::get_hdf5_type_ (type_t ());
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

  template <typename T>
  inline hdf5::private_::hdf5_buffer__
  hdf5_buffer (const T &t)
  {
    return hdf5::private_::hdf5_buffer__ (T (), &t, 1, 0, 0);
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
    : type (get_hdf5_type (T ()))
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
    hid_t                       type;
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
      static hid_t
      get_type (const T &h)
      {
        return hdf5_struct_helper <T::lhs_type::static_size>::get_type (h.lhs);
      }
      template <typename T>
      static hid_t
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
      static hid_t
      get_type (const T &h)
      {
        return h.type;
      }
      template <typename T>
      static hid_t
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
      : type (-1)
      , size (get_size <L>::static_size + get_size <R>::static_size)
      , data_ (0)
      {
        typedef hdf5_struct_holder <L, R> T;

        size += hdf5_struct_helper <T::static_size>::get_size (h);
        type  = get_hdf5_type (typename T::type_t ());

        typename T::type_t *buffer = new typename T::type_t [size];
        hdf5_struct_helper <T::static_size>::fill_buffer (h, buffer, 0);

        data_ = buffer;
      }

      void *
      data () const
      {
        return data_;
      }

      hid_t     type;
      size_t    size;
      void      *data_;
    };

  } // namespace private_
  } // namespace hdf5


  struct BS_API_PLUGIN hdf5_file;
  struct BS_API_PLUGIN hdf5_group_v2;
  struct BS_API_PLUGIN hdf5_storage_v2
  {
  private:

    static hdf5_storage_v2 *
    instance ();

    hdf5_storage_v2 ();
    ~hdf5_storage_v2 ();

    struct impl;
    impl *impl_;

    friend struct hdf5_file;
    friend struct hdf5_group_v2;
    friend struct impl;
  };

  struct BS_API_PLUGIN hdf5_file
  {
    hdf5_group_v2 
    operator [] (const std::string &name);

    hdf5_file (const std::string &file_name)
    : file_id_ (-1)
    , file_name_ (file_name)
    {
    }

    ~hdf5_file ();

  private:

    hid_t       file_id_;
    std::string file_name_;

    friend struct hdf5_group_v2;
    friend struct hdf5_property_v2;
    friend struct hdf5_storage_v2;
  };

  struct BS_API_PLUGIN hdf5_group_v2
  {
  private:
    typedef hdf5::private_::hdf5_buffer__ hdf5_buffer_t;
    typedef hdf5::private_::hdf5_pod__    hdf5_pod_t;
    typedef hdf5::private_::hdf5_struct__ hdf5_struct_t;

  public:

    hdf5_group_v2 (const hdf5_file &file, const std::string &name)
    : file_ (file)
    , name_ (name)
    {
    }

    template <typename L>
    hdf5_group_v2
    write (const char *dataset, const hdf5::private_::hdf5_value_holder_unary <L> &h)
    {
      return write_buffer (dataset, hdf5_buffer (h.v));
    }

    template <typename L, typename R>
    hdf5_group_v2 &
    write (const char *dataset, const hdf5::private_::hdf5_value_holder <L, R> &h)
    {
      using namespace hdf5::private_;
      return write_pod (dataset, hdf5_pod_impl <hdf5_value_holder <L, R> > (h));
    }

    hdf5_group_v2 &
    write (const char *dataset, const hdf5_buffer_t &buffer)
    {
      return write_buffer (dataset, buffer);
    }

    template <typename T>
    hdf5_group_v2 &
    write (const char *dataset, const shared_vector <T> &data)
    {
      return write_buffer (dataset, hdf5_buffer (data));
    }

    template <typename T>
    hdf5_group_v2 &
    write (const char *dataset, const hdf5_struct <T> &s)
    {
      return write_buffer (dataset, hdf5_buffer_t (s.data (), s.size, s.type));
    }

    template <typename L, typename R>
    hdf5_group_v2 &
    write (const char *dataset, const hdf5::private_::hdf5_struct_holder <L, R> &s)
    {
      return write_struct (dataset, hdf5_struct_t (s));
    }

  private:

    hdf5_group_v2 &
    write_buffer (const char *dataset, const hdf5_buffer_t &buffer);

    hdf5_group_v2 &
    write_pod (const char *dataset, const hdf5_pod_t &pod);

    hdf5_group_v2 &
    write_struct (const char *dataset, const hdf5_struct_t &s);

  private:
    hdf5_file   file_;
    std::string name_;
  };

  template <typename L, typename R>
  inline hdf5::private_::hdf5_value_holder <L, R>
  operator << (const hdf5::private_::hdf5_value_holder_unary <L> &lhs, const R &rhs)
  {
    return hdf5::private_::hdf5_value_holder <L, R> (lhs.v, rhs);
  }

  template <typename LL, typename LR, typename R>
  inline hdf5::private_::hdf5_value_holder <hdf5::private_::hdf5_value_holder <LL, LR>, R>
  operator << (const hdf5::private_::hdf5_value_holder <LL, LR> &lhs, const R &rhs)
  {
    return hdf5::private_::hdf5_value_holder <hdf5::private_::hdf5_value_holder <LL, LR>, R> (lhs, rhs);
  }

  template <typename T, typename R>
  inline hdf5::private_::hdf5_struct_holder <hdf5_struct <T>, R>
  operator << (const hdf5_struct <T> &lhs, const R &rhs)
  {
    return hdf5::private_::hdf5_struct_holder <hdf5_struct <T>, R> (lhs, rhs);
  }
  template <typename LL, typename LR, typename R>
  inline hdf5::private_::hdf5_struct_holder <hdf5::private_::hdf5_struct_holder <LL, LR>, R>
  operator << (const hdf5::private_::hdf5_struct_holder <LL, LR> &lhs, const R &rhs)
  {
    using namespace hdf5::private_;
    return hdf5_struct_holder <hdf5_struct_holder <LL, LR>, R> (lhs, rhs);
  }


} // namespace blue_sky

#endif // #ifndef BS_EAGLE_HDF5_STORAGE_V2_H_

