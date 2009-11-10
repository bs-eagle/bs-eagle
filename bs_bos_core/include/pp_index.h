/**
 * \file pp_index.h
 * \brief 
 * \author Sergey Miryanov
 * \date 15.10.2009
 * */
#ifndef BS_BOS_CORE_PP_INDEX_H_
#define BS_BOS_CORE_PP_INDEX_H_

namespace blue_sky {
namespace detail {

  template <size_t, bool, bool>
  struct pp_index
  {
  };

  template <bool is_w, bool is_g>
  struct pp_index <1, is_w, is_g>
  {
    enum
    {
      gas_po = 0,
      oil_po = 0,
      wat_po = 0,

      gas_sg = -1, gas_so = -1,
      oil_sg = -1, oil_so = -1,
      wat_sg = -1, wat_so = -1,
    };
  };

  template <>
  struct pp_index <2, true, false>
  {
    enum
    {
      oil_so = 0,
      oil_po = 1,
      wat_so = 2,
      wat_po = 3,

      gas_sg = -1, gas_so = -1, gas_po = -1,
      oil_sg = -1,
      wat_sg = -1,
    };
  };
  template <>
  struct pp_index <2, false, true>
  {
    enum
    {
      gas_sg = 0,
      gas_po = 1,
      oil_sg = 2,
      oil_po = 3,
      gas_so = 4,
      oil_so = 5,

      wat_sg = -1, wat_so = -1, wat_po = -1,
    };
  };

  template <bool is_w, bool is_g>
  struct pp_index <3, is_w, is_g>
  {
    enum
    {
      gas_sg = 0, gas_so = 1, gas_po = 2,
      oil_sg = 3, oil_so = 4, oil_po = 5,
      wat_sg = 6, wat_so = 7, wat_po = 8,
    };
  };

} // namespace detail
} // namespace blue_sky


#endif // #ifndef BS_BOS_CORE_PP_INDEX_H_
