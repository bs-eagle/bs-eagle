/**
 * \file auto_value.h
 * \brief helper to init scalar variables with auto value
 * \author Sergey Miryanov, idea Kodt@RSDN
 * \date 10.09.2008
 */
#ifndef BS_AUTO_VALUE_HELPER_H_
#define BS_AUTO_VALUE_HELPER_H_

namespace blue_sky
  {

  // TODO: type traits
  template<typename T, long value = 0>
  class auto_value
    {
    public:
      typedef T           data_t;
      typedef auto_value  &self_t;

      inline auto_value()
          : t_(T (value))
#ifdef _DEBUG
          , is_init_ (false)
#endif
      {
      }

      template< class V >
      inline auto_value(V v)
          : t_(v)
#ifdef _DEBUG
          , is_init_ (true)
#endif
      {
      }

      inline const T& data() const
        {
          return t_;
        }
      inline T& data()
      {
        return t_;
      }

      inline operator T  () const
        {
          return t_;
        }
      inline operator T& ()
      {
        return t_;
      }

#ifdef _DEBUG
      template< class V > inline self_t operator =   (V v)
      {
        t_ =   v;
        is_init_ = true;
        return *this;
      }
#else
      template< class V > inline self_t operator =   (V v)
      {
        t_ =   v;
        return *this;
      }
#endif
      template< class V > inline self_t operator +=  (V v)
      {
        t_ +=  v;
        return *this;
      }
      template< class V > inline self_t operator -=  (V v)
      {
        t_ -=  v;
        return *this;
      }
      template< class V > inline self_t operator *=  (V v)
      {
        t_ *=  v;
        return *this;
      }
      template< class V > inline self_t operator /=  (V v)
      {
        t_ /=  v;
        return *this;
      }
      template< class V > inline self_t operator %=  (V v)
      {
        t_ %=  v;
        return *this;
      }
      template< class V > inline self_t operator &=  (V v)
      {
        t_ &=  v;
        return *this;
      }
      template< class V > inline self_t operator |=  (V v)
      {
        t_ |=  v;
        return *this;
      }
      template< class V > inline self_t operator ^=  (V v)
      {
        t_ ^=  v;
        return *this;
      }
      template< class V > inline self_t operator <<= (V v)
      {
        t_ <<= v;
        return *this;
      }
      template< class V > inline self_t operator >>= (V v)
      {
        t_ >>= v;
        return *this;
      }

      inline self_t operator ++ ()
      {
        ++t_;
        return *this;
      }
      inline data_t operator ++ (int)
      {
        return t_++;
      }
      inline self_t operator -- ()
      {
        --t_;
        return *this;
      }
      inline data_t operator -- (int)
      {
        return t_--;
      }

#ifdef _DEBUG
      inline bool is_init () const
        {
          return is_init_;
        }
#endif

    private:
      T         t_;
#ifdef _DEBUG
      bool      is_init_;
#endif
    };

  //template<>
  //class auto_value <float, 0>
  //{
  //public:
  //  typedef float       T;
  //  typedef T           data_t;
  //  typedef auto_value  &self_t;

  //  enum { value = 0 };

  //  inline auto_value()
  //    : t_(value)
  //  {
  //  }

  //  template< class V >
  //  inline auto_value(V v)
  //    : t_(v)
  //  {
  //  }

  //  inline const T& data() const  { return t_; }
  //  inline T& data()              { return t_; }

  //  inline operator T  () const   { return t_; }
  //  inline operator T& ()         { return t_; }

  //  template< class V > inline self_t operator =   (V v)  { t_ =   v; return *this; }
  //  template< class V > inline self_t operator +=  (V v)  { t_ +=  v; return *this; }
  //  template< class V > inline self_t operator -=  (V v)  { t_ -=  v; return *this; }
  //  template< class V > inline self_t operator *=  (V v)  { t_ *=  v; return *this; }
  //  template< class V > inline self_t operator /=  (V v)  { t_ /=  v; return *this; }
  //  template< class V > inline self_t operator %=  (V v)  { t_ %=  v; return *this; }
  //  template< class V > inline self_t operator &=  (V v)  { t_ &=  v; return *this; }
  //  template< class V > inline self_t operator |=  (V v)  { t_ |=  v; return *this; }
  //  template< class V > inline self_t operator ^=  (V v)  { t_ ^=  v; return *this; }
  //  template< class V > inline self_t operator <<= (V v)  { t_ <<= v; return *this; }
  //  template< class V > inline self_t operator >>= (V v)  { t_ >>= v; return *this; }

  //  inline self_t operator ++ ()                          { ++t_; return *this;     }
  //  inline data_t operator ++ (int)                       { return t_++;            }
  //  inline self_t operator -- ()                          { --t_; return *this;     }
  //  inline data_t operator -- (int)                       { return t_--;            }

  //private:
  //  T t_;
  //};

  //template<>
  //class auto_value <double>
  //{
  //public:
  //  typedef float       T;
  //  typedef T           data_t;
  //  typedef auto_value  &self_t;

  //  enum { value = 0 };

  //  inline auto_value()
  //    : t_(value)
  //  {
  //  }

  //  template< class V >
  //  inline auto_value(V v)
  //    : t_(v)
  //  {
  //  }

  //  inline const T& data() const  { return t_; }
  //  inline T& data()              { return t_; }

  //  inline operator T  () const   { return t_; }
  //  inline operator T& ()         { return t_; }

  //  template< class V > inline self_t operator =   (V v)  { t_ =   v; return *this; }
  //  template< class V > inline self_t operator +=  (V v)  { t_ +=  v; return *this; }
  //  template< class V > inline self_t operator -=  (V v)  { t_ -=  v; return *this; }
  //  template< class V > inline self_t operator *=  (V v)  { t_ *=  v; return *this; }
  //  template< class V > inline self_t operator /=  (V v)  { t_ /=  v; return *this; }
  //  template< class V > inline self_t operator %=  (V v)  { t_ %=  v; return *this; }
  //  template< class V > inline self_t operator &=  (V v)  { t_ &=  v; return *this; }
  //  template< class V > inline self_t operator |=  (V v)  { t_ |=  v; return *this; }
  //  template< class V > inline self_t operator ^=  (V v)  { t_ ^=  v; return *this; }
  //  template< class V > inline self_t operator <<= (V v)  { t_ <<= v; return *this; }
  //  template< class V > inline self_t operator >>= (V v)  { t_ >>= v; return *this; }

  //  inline self_t operator ++ ()                          { ++t_; return *this;     }
  //  inline data_t operator ++ (int)                       { return t_++;            }
  //  inline self_t operator -- ()                          { --t_; return *this;     }
  //  inline data_t operator -- (int)                       { return t_--;            }

  //private:
  //  T t_;
  //};

} // namespace blue_sky


#endif  // #ifndef BS_AUTO_VALUE_HELPER_H_
