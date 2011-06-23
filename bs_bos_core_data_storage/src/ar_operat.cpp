#include "bs_bos_core_data_storage_stdafx.h"
#include "ar_operat.h"

namespace blue_sky
  {
  ar_operat::ar_operat()
  {
    init();
  }

  ar_operat::ar_operat(const std::string &tname)
  {
    init ();
    name = tname;
  }

  ar_operat::ar_operat(operat_func toper, int tpriority, int tnum_of_arg, const std::string &tname)
      : oper(toper), priority(tpriority), num_of_arg(tnum_of_arg), name(tname)
  {}

  ar_operat::~ar_operat()
  {
    init();
  }

  void ar_operat::init()
  {
    oper = 0;
    priority = 0;
    num_of_arg = 0;
    name = "";
  }

  int ar_operat::set_name (const std::string &tname)
  {
    if (!tname.length())
      {
        BS_ASSERT (tname != "");
        throw bs_exception ("ar_operat::set_name", "Invalid argument");
      }
    name = tname;
    return 0;
  }

  const std::string &
  ar_operat::get_name() const
  {
    return name;
  }

  int ar_operat::operator==(const std::string &tname) const
  {
    return strcmp(name.c_str(),tname.c_str());
  }

  bool ar_operat::operator<(const ar_operat &op) const
  {
    return (name < op.name);
  }
}
