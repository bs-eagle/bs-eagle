/**
 * \file py_read_class.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 28.10.2009
 * */
#include "bs_bos_core_data_storage_stdafx.h"
#include "read_class.h"
#include "py_read_class.h"

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky {
namespace python {

  template <typename T>
  std::string
  read_line (T *t)
  {
    char buf[CHAR_BUF_LEN + 1] = {0};
    int len = t->read_line (buf, CHAR_BUF_LEN);
    if (len <= 0)
      {
        bs_throw_exception (boost::format ("Can't read line from file %s") % t->get_prefix ());
      }

    return std::string (buf, len);
  }

  PY_EXPORTER (fread_exporter, default_exporter)
    .def ("read_line", read_line <T>)
  PY_EXPORTER_END;

  void
  export_FRead ()
  {
    base_exporter <FRead, fread_exporter>::export_class ("data_reader");
  }


} // namespace python
} // namespace blue_sky
#endif