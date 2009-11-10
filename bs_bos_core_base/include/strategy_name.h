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
  struct strategy_name <base_strategy_di>
  {
    static std::string
    name () 
    {
      return "di";
    }
  };

  template <>
  struct strategy_name <base_strategy_fi>
  {
    static std::string
    name () 
    {
      return "fi";
    }
  };

  template <>
  struct strategy_name <base_strategy_mixi>
  {
    static std::string
    name () 
    {
      return "mixi";
    }
  };

}   // namespace tools
}   // namespace blue_sky


#endif  // #ifndef BS_BOS_CORE_BASE_STRATEGY_NAME_H_
