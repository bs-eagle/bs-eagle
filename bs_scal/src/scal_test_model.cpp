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

    test_model::test_model (bs_type_ctor_param param /* = NULL */)
        : n_phases (0)
        , cell_count (0)
    {

    }

    test_model::test_model (const this_t& tm)
    : bs_refcounter (tm), objbase (tm)
    {
      BS_ASSERT (false && "NOT IMPL YET");
    }

    void
    test_model::init (int np, int cc)
    {
      n_phases = np;
      cell_count = cc;

      relative_perm.assign (cell_count * n_phases, 0);
      s_deriv_relativ_perm.assign (cell_count * 2 * (n_phases - 1), 0);

      cap_pressure.assign (cell_count * (n_phases - 1), 0);
      s_deriv_cap_pressure.assign (cell_count * (n_phases - 1), 0);

      sat.assign (cell_count * (n_phases - 1), 0);
    }

    test_model::item_t *
    test_model::get_relative_perm (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &relative_perm[cell_index * n_phases];
    }

    
    test_model::item_t *
    test_model::get_s_deriv_relativ_perm (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &s_deriv_relativ_perm[cell_index * 2 * (n_phases - 1)];
    }

    test_model::item_t *
    test_model::get_cap_pressure (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &cap_pressure[cell_index * (n_phases - 1)];
    }

    test_model::item_t *
    test_model::get_s_deriv_cap_pressure (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &s_deriv_cap_pressure[cell_index * (n_phases - 1)];
    }

    test_model::item_t *
    test_model::get_sat (int cell_index)
    {
      BS_ASSERT (cell_index >= 0 && cell_index < cell_count);

      return &sat [cell_index * (n_phases - 1)];
    }

    int
    test_model::get_cell_count () const
      {
        return cell_count;
      }

    //////////////////////////////////////////////////////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (test_model);
  BLUE_SKY_TYPE_STD_COPY (test_model);

  BLUE_SKY_TYPE_IMPL(test_model,  objbase, "test_model", "test_model SCAL calculation class", "test_model SCAL calculation");

  } // namespace scal
} // namespace blue_sky
