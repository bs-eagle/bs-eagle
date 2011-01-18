/** 
 * \file strategy_name.h
 * \brief helper to obtain name of strategy
 * \author Sergey Miryanov
 * \date 14.07.2009
 * */
#ifndef BS_BOS_CORE_BASE_STRATEGY_NAME_H_
#define BS_BOS_CORE_BASE_STRATEGY_NAME_H_

namespace blue_sky {
namespace tools {

  template <typename strategy_t>
  struct strategy_name
  {
  };

  template <>
  struct strategy_name <base_strategy_did>
  {
    static std::string
    name () 
    {
      return "did";
    }
  };

  template <>
  struct strategy_name <base_strategy_fif>
  {
    static std::string
    name () 
    {
      return "fif";
    }
  };

  template <>
  struct strategy_name <base_strategy_dif>
  {
    static std::string
    name () 
    {
      return "dif";
    }
  };
  
  
  template <>
  struct strategy_name <base_strategy_dld>
  {
    static std::string
    name () 
    {
      return "dld";
    }
  };

  template <>
  struct strategy_name <base_strategy_flf>
  {
    static std::string
    name () 
    {
      return "flf";
    }
  };

  template <>
  struct strategy_name <base_strategy_dlf>
  {
    static std::string
    name () 
    {
      return "dlf";
    }
  };

}   // namespace tools
}   // namespace blue_sky


#endif  // #ifndef BS_BOS_CORE_BASE_STRATEGY_NAME_H_
