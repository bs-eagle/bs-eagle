/**
 * \file well_rate_call_proxy.h
 * \brief helper class to call functors in connection loop
 * \author Sergey Miryanov
 * \date 21.11.2008
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_CALL_PROXY_H_
#define BS_WELLS_WELL_RATE_CONTROL_CALL_PROXY_H_

namespace blue_sky
  {

  template <typename wrapped_t, typename callee_t>
  struct one_call_proxy
    {
      one_call_proxy (wrapped_t *wrapped, callee_t callee)
          : wrapped_ (wrapped)
          , callee_ (callee)
      {
      }

      template <typename locked_connection_t, typename data_t, typename params_t>
      void
      operator () (const locked_connection_t &c, const data_t &data, params_t &params) const
        {
          ((*wrapped_).*callee_) (c, data, params);
        }
      template <typename locked_connection_t, typename data_t>
      void
      operator () (const locked_connection_t &c, const data_t &data) const
        {
          ((*wrapped_).*callee_) (c, data);
        }
      template <typename locked_connection_t, typename params_t>
      void
      operator () (const locked_connection_t &c, params_t &params) const
        {
          ((*wrapped_).*callee_) (c, params);
        }

      wrapped_t *wrapped_;
      callee_t  callee_;
    };

  template <typename wrapped_t, typename callee_t>
  struct two_call_proxy
    {
      two_call_proxy (wrapped_t *wrapped, callee_t first_callee, callee_t second_callee)
          : wrapped_ (wrapped)
          , first_callee_ (first_callee)
          , second_callee_ (second_callee)
      {
      }

      template <typename locked_connection_t, typename data_t, typename params_t>
      void
      operator () (const locked_connection_t &c, const data_t &data, params_t &params) const
        {
          ((*wrapped_).*first_callee_) (c, data, params);
          ((*wrapped_).*second_callee_) (c, data, params);
        }

      wrapped_t *wrapped_;
      callee_t  first_callee_;
      callee_t  second_callee_;
    };

  template <typename this_t, typename callee_t>
  static one_call_proxy <this_t, callee_t>
  one_call (this_t *this_, callee_t callee)
  {
    return one_call_proxy <this_t, callee_t> (this_, callee);
  }
  template <typename this_t, typename callee_t>
  static two_call_proxy <this_t, callee_t>
  two_call (this_t *this_, callee_t first_callee, callee_t second_callee)
  {
    return two_call_proxy <this_t, callee_t> (this_, first_callee, second_callee);
  }

} // namespace blue_sky

#endif  // #ifndef BS_WELL_WELL_RATE_CONTROL_CALL_PROXY_H_

