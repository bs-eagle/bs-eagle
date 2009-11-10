/**
 * \file scal_test_model.cpp
 * \brief impl of scal::test_model
 * \author Sergey Miryanov
 * \date 23.05.2008
 * */
#include "bs_scal_stdafx.h"
#include "scal_test_model.h"

namespace blue_sky
  {
  namespace scal
    {

    template <typename strategy_t>
    test_model<strategy_t>::test_model (bs_type_ctor_param param /* = NULL */)
        : n_phases (0)
        , cell_count (0)
    {

    }

    template <typename strategy_t>
    test_model<strategy_t>::test_model (const this_t& tm)
    : bs_refcounter (tm), objbase (tm)
    {
      BS_ASSERT (false && "NOT IMPL YET");
    }

    template <typename strategy_t>
    void
    test_model<strategy_t>::init (int np, int cc)
    {
      n_phases = np;
      cell_count = cc;

      relative_perm.assign (cell_count * n_phases, 0);
      s_deriv_relativ_perm.assign (cell_count * 2 * (n_phases - 1), 0);

      cap_pressure.assign (cell_count * (n_phases - 1), 0);
      s_deriv_cap_pressure.assign (cell_count * (n_phases - 1), 0);

      sat.assign (cell_count * (n_phases - 1), 0);
    }

    template <typename strategy_t>
    typename test_model<strategy_t>::item_t *
    test_model<strategy_t>::get_relative_perm (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &relative_perm[cell_index * n_phases];
    }

    template <typename strategy_t>
    typename test_model<strategy_t>::item_t *
    test_model<strategy_t>::get_s_deriv_relativ_perm (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &s_deriv_relativ_perm[cell_index * 2 * (n_phases - 1)];
    }

    template <typename strategy_t>
    typename test_model<strategy_t>::item_t *
    test_model<strategy_t>::get_cap_pressure (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &cap_pressure[cell_index * (n_phases - 1)];
    }

    template <typename strategy_t>
    typename test_model<strategy_t>::item_t *
    test_model<strategy_t>::get_s_deriv_cap_pressure (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &s_deriv_cap_pressure[cell_index * (n_phases - 1)];
    }

    template <typename strategy_t>
    typename test_model<strategy_t>::item_t *
    test_model<strategy_t>::get_sat (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &sat [cell_index * (n_phases - 1)];
    }

    template <typename strategy_t>
    int
    test_model<strategy_t>::get_cell_count () const
      {
        return cell_count;
      }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE_T_DEF(test_model,(class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF(test_model,(class));
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (test_model<base_strategy_fi>), 1, (objbase), "test_model_fi", "test_model_fi", "test_model_fi", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (test_model<base_strategy_di>), 1, (objbase), "test_model_di", "test_model_di", "test_model_di", false);

  } // namespace scal
} // namespace blue_sky
