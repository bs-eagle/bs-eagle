/**
 *       \file  facility_base.h
 *      \brief  Base class (interface) for facilities (wells, unnamed wells e.t.c)
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Should be renamed to facility_iface.h
 * */
#ifndef BS_FACILITY_BASE_H_
#define BS_FACILITY_BASE_H_

#include "calc_model.h"

namespace blue_sky {

  /**
   * \class facility_base
   * \brief Base class (interface) for facilities (wells e.t.c)
   * \todo  Should be renamed to facility_iface
   * */
  class BS_API_PLUGIN facility_base : public objbase
  {
  public:
    typedef v_long                          index_array_t;
    typedef v_double                        item_array_t;
    typedef v_double                        rhs_item_array_t;
    typedef t_long                          index_t;

    typedef rs_mesh_iface                   mesh_iface_t;
    typedef jac_matrix_iface                jmatrix_t;

    typedef smart_ptr <calc_model, true>    sp_calc_model_t;
    typedef smart_ptr <mesh_iface_t, true>  sp_mesh_iface_t;
    typedef smart_ptr <jmatrix_t, true>     sp_jmatrix_t;

  public:
    //! destructor
    virtual ~facility_base () {}

    /**
     * \brief  Fills Jacobian rows
     * \param  rows Array of Jacobian Rows
     * */
    virtual void 
    fill_rows (index_array_t &rows) const = 0;

    /**
     * \brief  Fills Jacobian colls and values, uses eliminate for 
     *         fill values
     * \param  dt
     * \param  block_size
     * \param  rows
     * \param  cols
     * \param  values
     * \param  markers
     * */
    virtual void 
    fill_jacobian (double dt, index_t block_size, const index_array_t &rows, index_array_t &cols, rhs_item_array_t &values, index_array_t &markers) const = 0;

    /**
     * \brief  Fills rhs array with rate values
     * \param  dt
     * \param  n_phases
     * \param  is_g
     * \param  is_o
     * \param  is_w
     * \param  rhs Array of rhs values
     * */
    virtual void 
    fill_rhs (double dt, index_t n_phases, bool is_g, bool is_o, bool is_w, rhs_item_array_t &rhs) const = 0;

    /**
     * \brief  Calculates rate and deriv values for well and 
     *         well perforations (connections)
     * \param  is_start 
     * \param  dt
     * \param  calc_model
     * \param  mesh
     * \param  jmatrix
     * */
    virtual void 
    process (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix) = 0;

    /**
     * \brief  Restores solution
     * \param  dt
     * \param  p_sol primary solution vector
     * \param  s_sol secondary solution vector
     * \param  block_size size of one block in vectors
     * */
    virtual void 
    restore_solution (double dt, const item_array_t &p_sol, const item_array_t &s_sol, index_t block_size) = 0;

    /**
     * \brief  Performs actions before start of each large step
     * \param  calc_model
     * \param  mesh
     * */
    virtual void 
    pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) = 0;

    /**
     * \brief  Performs actions before start of each small step
     * */
    virtual void 
    pre_small_step () = 0;

    /**
     * \brief  Performs actions before start of each newton step
     * */
    virtual void 
    pre_newton_step () = 0;

    /**
     * \brief  Restores 'internal-state' of well if small step failed
     * */
    virtual void 
    restart_small_step () = 0;

    /**
     * \brief  Restores 'internal-state' of well if newton step failed
     * */
    virtual void 
    restart_newton_step () = 0;
  };

} // namespace blue_sky


#endif  // #ifndef BS_FACILITY_BASE_H_

