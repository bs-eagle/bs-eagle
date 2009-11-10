/**
 * \file default_well.h
 * \brief default well for bs_bos_core
 * \author Sergey Miryanov
 * \date 20.05.2009
 * */
#ifndef BS_BOS_CORE_DEFAULT_WELL_H_
#define BS_BOS_CORE_DEFAULT_WELL_H_

#include "calc_well.h"
#include "default_connection.h"

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  class BS_API_PLUGIN default_well : public well <strategy_t>
  {
  public:
      typedef well <strategy_t>                       	base_t;
      typedef well <strategy_t>                         well_t;
      typedef typename base_t::item_array_t					    item_array_t;
      typedef typename base_t::rhs_item_array_t			    rhs_item_array_t;
      typedef typename base_t::index_array_t				    index_array_t;
      typedef typename base_t::index_t							    index_t;
      typedef typename base_t::item_t								    item_t;
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
      typedef std::map <connection_block_t, connection_index_t>   connection_map_t;

  public:

    default_well (const std::string &well_name);

    BLUE_SKY_TYPE_DECL_T (default_well <strategy_t>);

    void
    restore_solution (double dt, const item_array_t &p_sol, const item_array_t &s_sol, index_t block_size);

    void
    process (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

    void
    eliminate (rhs_item_t *array, index_t rw_index, index_t wr_index, double dt, index_t block_size) const;

    void
    clear_data ();

    array_ext <item_t>
    get_ww_value ();

    array_ext <item_t>
    get_bw_value ();

  protected:
    void
    process_internal (bool is_start, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

    template <bool is_prod>
    void
    calc_rate_and_derivs (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

    template <bool is_w, bool is_g, bool is_o, bool is_prod>
    void
    calc_rate_and_derivs_concrete (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

  public:
    sp_connection_t
    add_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block);

    sp_connection_t
    get_connection (index_t idx) const;

    sp_connection_t
    get_connection_map (index_t n_block) const;

    virtual size_t 
    get_connections_count () const;

  public:

    item_t                      ww_value;
    item_t                      bw_value;

    connection_list_t           connection_list_;
    connection_map_t            connection_map_;
  };

  bool
  default_well_register_types (const blue_sky::plugin_descriptor &pd);


} // namespace wells
} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_DEFAULT_WELL_H_

