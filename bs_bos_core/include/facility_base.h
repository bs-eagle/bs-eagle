/**
 * \file facility_base.h
 * \brief base class for facilities (wells and unnamed wells)
 * \author Sergey Miryanov
 * \date 16.07.2008
 */
#ifndef BS_FACILITY_BASE_H_
#define BS_FACILITY_BASE_H_

#include "calc_model.h"

namespace blue_sky {

  template <typename strategy_t>
  class BS_API_PLUGIN facility_base : public objbase
  {
  public:
    typedef typename strategy_t::index_array_t    index_array_t;
    typedef typename strategy_t::item_array_t     item_array_t;
    typedef typename strategy_t::rhs_item_array_t rhs_item_array_t;
    typedef typename strategy_t::index_t          index_t;

    typedef calc_model <strategy_t>               calc_model_t;
    typedef rs_mesh_iface <strategy_t>            mesh_iface_t;
    typedef jacobian_matrix <strategy_t>          jmatrix_t;

    typedef smart_ptr <calc_model_t, true>        sp_calc_model_t;
    typedef smart_ptr <mesh_iface_t, true>        sp_mesh_iface_t;
    typedef smart_ptr <jmatrix_t, true>           sp_jmatrix_t;

  public:
    //! destructor
    virtual ~facility_base () {}

    virtual void fill_rows (index_array_t &rows) const = 0;
    virtual void fill_jacobian (double dt, index_t block_size, const index_array_t &rows, index_array_t &cols, rhs_item_array_t &values, index_array_t &markers) const = 0;
    virtual void fill_rhs (double dt, index_t n_phases, bool is_g, bool is_o, bool is_w, rhs_item_array_t &rhs) const = 0;

    virtual void process (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix) = 0;

    virtual void restore_solution (double dt, const item_array_t &p_sol, const item_array_t &s_sol, index_t block_size) = 0;

    virtual void pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) = 0;
    virtual void pre_small_step () = 0;
    virtual void pre_newton_step () = 0;

    virtual void restart_small_step () = 0;
    virtual void restart_newton_step () = 0;
  };

} // namespace blue_sky


#endif  // #ifndef BS_FACILITY_BASE_H_

