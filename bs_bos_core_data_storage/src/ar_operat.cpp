#include "bs_bos_core_data_storage_stdafx.h"
#include "ar_operat.h"

namespace blue_sky
  {
	template <typename strategy_t>
  ar_operat<strategy_t>::ar_operat()
  {
    init();
  }

	template <typename strategy_t>
  ar_operat<strategy_t>::ar_operat(const std::string &tname)
  {
    init ();
    name = tname;
  }

	template <typename strategy_t>
  ar_operat<strategy_t>::ar_operat(operat_func toper, int tpriority, int tnum_of_arg, const std::string &tname)
      : oper(toper), priority(tpriority), num_of_arg(tnum_of_arg), name(tname)
  {}

	template <typename strategy_t>
  ar_operat<strategy_t>::~ar_operat()
  {
    init();
  }

	template <typename strategy_t>
  void ar_operat<strategy_t>::init()
  {
    oper = 0;
    priority = 0;
    num_of_arg = 0;
    name = "";
  }

	template <typename strategy_t>
  int ar_operat<strategy_t>::set_name (const std::string &tname)
  {
    if (!tname.length())
      {
        BS_ASSERT (tname != "");
        throw bs_exception ("ar_operat::set_name", "Invalid argument");
      }
    name = tname;
    return 0;
  }

	template <typename strategy_t>
  const std::string &
  ar_operat<strategy_t>::get_name() const
  {
    return name;
  }

	template <typename strategy_t>
  int ar_operat<strategy_t>::operator==(const std::string &tname) const
  {
    return strcmp(name.c_str(),tname.c_str());
  }

	template <typename strategy_t>
  bool ar_operat<strategy_t>::operator<(const ar_operat &op) const
  {
    return (name < op.name);
  }

	template class ar_operat <base_strategy_fi>;
  template class ar_operat <base_strategy_di>;
	template class ar_operat <base_strategy_mixi>;
}
