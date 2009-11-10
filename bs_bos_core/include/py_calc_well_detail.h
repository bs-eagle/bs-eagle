/**
 * \file py_calc_well_detail.h
 * \brief
 * \author Sergey Miryanov
 * \date 01.07.2009
 * */
#ifndef BS_BOS_CORE_PY_CALC_WELL_DETAIL_H_
#define BS_BOS_CORE_PY_CALC_WELL_DETAIL_H_

namespace blue_sky {
namespace python {
namespace detail {

  template <typename connection_t>
  static void
  set_connection_density (connection_t *c, double density)
  {
    c->density = density;
  }

  template <typename well_t>
  static typename well_t::index_t
  get_well_i_coord (const well_t *well)
  {
    return well->i_coord_;
  }
  template <typename well_t>
  static typename well_t::index_t
  get_well_j_coord (const well_t *well)
  {
    return well->j_coord_;
  }
  template <typename well_t>
  static void
  set_well_i_coord (well_t *well, typename well_t::index_t i)
  {
    well->i_coord_ = i;
  }
  template <typename well_t>
  static void
  set_well_j_coord (well_t *well, typename well_t::index_t j)
  {
    well->j_coord_ = j;
  }

  template <typename well_t>
  static typename well_t::item_t
  get_well_wefac (const well_t *well)
  {
    return well->exploitation_factor_;
  }

  template <typename well_t>
  static typename well_t::index_t
  get_well_state (const well_t *well)
  {
    return well->well_state_.state;
  }

  template <typename well_t>
  static bool
  get_well_is_work (const well_t *well)
  {
    return well->well_state_.is_work;
  }
  template <typename well_t>
  static void
  set_well_is_work (well_t *well, bool is_work)
  {
    well->well_state_.is_work = is_work;
  }
  template <typename well_t>
  static bool
  well_is_production (const well_t *well)
  {
    return well->get_controller ()->is_production ();
  }

} // namespace detail
} // namespace python
} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_PY_CALC_WELL_DETAIL_H_



