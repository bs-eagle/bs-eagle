/**
 * \file named_pbase_type_helper.h
 * \brief 
 * \author Sergey Miryanov
 * \date 15.07.2009
 * */
#ifndef BS_BOS_CORE_BASE_NAMED_PBASE_TYPE_HELPER_H_
#define BS_BOS_CORE_BASE_NAMED_PBASE_TYPE_HELPER_H_

namespace blue_sky {
namespace tools {

  template <size_t>
  struct named_pbase_value_type_helper
  {
  };

  template <>
  struct named_pbase_value_type_helper <named_pbase::PT_INT>
  {
    typedef int type;

    template <typename T, typename idx_t>
    static type
    get (const T *t, idx_t idx)
    {
      return t->get_int (idx);
    }

    template <typename T, typename idx_t>
    static type
    get_d (const T *t, idx_t idx, const type &def)
    {
      return t->get_int_d (idx, def);
    }
  };

  template <>
  struct named_pbase_value_type_helper <named_pbase::PT_FLOAT>
  {
    typedef double type;

    template <typename T, typename idx_t>
    static type
    get (const T *t, idx_t idx)
    {
      return t->get_float (idx);
    }

    template <typename T, typename idx_t>
    static type
    get_d (const T *t, idx_t idx, const type &def)
    {
      return t->get_float_d (idx, def);
    }
  };

  template <>
  struct named_pbase_value_type_helper <named_pbase::PT_STR>
  {
    typedef std::string type;

    template <typename T, typename idx_t>
    static type
    get (const T *t, idx_t idx)
    {
      return t->get_str (idx);
    }

    template <typename T, typename idx_t>
    static type
    get_d (const T *t, idx_t idx, const type &def)
    {
      return t->get_str_d (idx, def);
    }
  };

} // namespace tools
} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_BASE_NAMED_PBASE_TYPE_HELPER_H_