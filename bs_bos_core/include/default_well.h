/**
 *       \file  default_well.h
 *      \brief  Default implementation of well
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  20.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Should be moved to src/
 * */
#ifndef BS_BOS_CORE_DEFAULT_WELL_H_
#define BS_BOS_CORE_DEFAULT_WELL_H_

#include "calc_well.h"
#include "default_connection.h"

namespace blue_sky {
namespace wells {

  /**
   * \class default_well
   * \brief Default implementation of well
   * */
  template <typename strategy_t>
  class BS_API_PLUGIN default_well : public well <strategy_t>
  {
  public:
      typedef well <strategy_t>                         base_t;
      typedef well <strategy_t>                         well_t;
      typedef typename base_t::item_array_t             item_array_t;
      typedef typename base_t::rhs_item_array_t         rhs_item_array_t;
      typedef typename base_t::index_array_t            index_array_t;
      typedef typename base_t::index_t                  index_t;
      typedef typename base_t::item_t                   item_t;
      typedef typename base_t::rhs_item_t               rhs_item_t;

      typedef typename base_t::sp_calc_model_t          sp_calc_model_t;
      typedef typename base_t::sp_mesh_iface_t          sp_mesh_iface_t;
      typedef typename base_t::sp_jmatrix_t             sp_jmatrix_t;

      typedef typename base_t::connection_t             connection_t;
      typedef default_connection <strategy_t>           default_connection_t;

      typedef typename base_t::sp_connection_t          sp_connection_t;
      typedef smart_ptr <default_connection_t>          sp_default_connection_t;

      typedef seq_vector <sp_default_connection_t>      connection_list_t;
      typedef index_t                                   connection_block_t;
      typedef index_t                                   connection_index_t;

      typedef std::map <connection_block_t, connection_index_t>   connection_map_t;   //!< \todo Obsolete, deprecated

  public:

    /**
     * \brief  ctor
     * \todo   Obsolete, deprecated
     * \param  well_name Name of well
     * */
    default_well (const std::string &well_name);

    //! blue-sky type declaration
    BLUE_SKY_TYPE_DECL_T (default_well <strategy_t>);

    /**
     * \brief  Restores solution
     * \param  dt
     * \param  p_sol primary solution vector
     * \param  s_sol secondary solution vector
     * \param  block_size size of one block in vectors
     * */
    void
    restore_solution (double dt, const item_array_t &p_sol, const item_array_t &s_sol, index_t block_size);

    /**
     * \brief  Calculates rate and deriv values for well and 
     *         well perforations (connections)
     * \param  is_start 
     * \param  dt
     * \param  calc_model
     * \param  mesh
     * \param  jmatrix
     * */
    void
    process (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

    /**
     * \brief  Calculates Jacobian value and stores it in array
     * \param  array Array of Jacobian values
     * \param  rw_index
     * \param  wr_index
     * \param  dt
     * \param  block_size
     * */
    void
    eliminate (rhs_item_t *array, index_t rw_index, index_t wr_index, double dt, index_t block_size) const;

    /**
     * \brief  Clears well and well perforations data
     * */
    void
    clear_data ();

    /**
     * \brief  Returns ww_value
     * \todo   Obsolete, should be removed
     * */
    array_ext <item_t>
    get_ww_value ();

    /**
     * \brief  Returns bw_value
     * \todo   Obsolete, should be removed 
     * */
    array_ext <item_t>
    get_bw_value ();

  protected:
    /**
     * \brief  Calculates rate and deriv values for well and 
     *         well perforations (connections)
     * \param  is_start 
     * \param  dt
     * \param  calc_model
     * \param  mesh
     * \param  jmatrix
     * */
    void
    process_internal (bool is_start, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

    /**
     * \brief  Calculates rate and derivs for well perforations,
     *         parametrized with is_prod value. is_prod true for
     *         production wells.
     * \param  calc_model
     * \param  mesh
     * \param  jmatrix
     * */
    template <bool is_prod>
    void
    calc_rate_and_derivs (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

    /**
     * \brief  Calculates rate and derivs for well perforations,
     *         parametrized with is_w, is_g, is_o (for water phase, 
     *         gas phase, oil phase) and is_prod (for production 
     *         well). Called from calc_rate_and_derivs.
     * \param  calc_model
     * \param  mesh
     * \param  jmatrix
     * */
    template <bool is_w, bool is_g, bool is_o, bool is_prod>
    void
    calc_rate_and_derivs_concrete (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

  public:
    /**
     * \brief  Adds connection (perforation) to well and return it
     * \param  i_coord i coordinate of perforation
     * \param  j_coord j coordinate of perforation
     * \param  k_coord k coordinate of perforation
     * \param  n_block Index of block (cell) in mesh for (i, j, k) coordinates
     * \return Created connection
     * */
    sp_connection_t
    add_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block);

    /**
     * \brief  Returns connection (perforation) with index idx
     * \param  idx Index of connection
     * \return connection instance on success otherwise null pointer
     * */
    sp_connection_t
    get_connection (index_t idx) const;

    /**
     * \brief  Returns connection (perforation) with n_block
     * \param  n_block Value of block index
     * \return connection instance on success otherwise null pointer
     * */
    sp_connection_t
    get_connection_map (index_t n_block) const;

    /**
     * \brief  Returns count of connection (perforation) for well
     * \return Count of connections
     * */
    virtual size_t 
    get_connections_count () const;

  public:

    item_t                      ww_value;           //!< WW value
    item_t                      bw_value;           //!< BW value

    connection_list_t           connection_list_;   //!< List of connections (default_connection)
    connection_map_t            connection_map_;    //!< \todo Obsolete
  };

  /**
   * \brief  Registers default_well and default_connection
   *         types in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return True if all types registered successfully
   * */
  bool
  default_well_register_types (const blue_sky::plugin_descriptor &pd);


} // namespace wells
} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_DEFAULT_WELL_H_

