/**
 * \file str_functor.h
 * \brief functor that used in event_manager and event_base parse functions. it should be used instead std or boost function to avoid ICE under MSVS
 * \author Morozov Andrey
 * \date 06.08.2008
 * */
#ifndef BS_STR_FUNCTOR_H_
#define BS_STR_FUNCTOR_H_

namespace blue_sky
  {

  //!functor for boost spirit string parsers handling
  template <typename T>
  struct str_functor
    {
      typedef void (T::*method_t)(const char*, const char*);
      typedef T self_t;

      str_functor (self_t *self, method_t method)
          : self (self)
          , method (method)
      {
      }

      void operator()(const char *first, const char *last) const
        {
          (self->*method)(first, last);
        }

      self_t *self;
      method_t method;
    };

} // namespace blue_sky

#endif  // #ifndef BS_STR_FUNCTOR_H_

