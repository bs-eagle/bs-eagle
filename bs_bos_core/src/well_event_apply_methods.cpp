/**
 * \file well_event_apply_methods.cpp
 * \brief impl of <well_event>::apply methods
 * \author Sergey Miryanov
 * \date 11.07.2008
 * */
#include "stdafx.h"

#include "reservoir.h"
#include "reservoir_simulator.h"
#include "well_events.h"
#include "calc_well.h"
#include "well_controller.h"
#include "well_limit_operation.h"
#include "well_connection.h"
#include "facility_manager.h"

#include "calc_model.h"

namespace blue_sky
  {


#define DECL_GET_WELL_EVENT_NAMES(event_name)                           \
  template <typename strategy_t>                                        \
  BS_API_PLUGIN std::string                                             \
  event_name <strategy_t>::get_well_name () const                       \
  {                                                                     \
  return main_params_->get_WELL_NAME ("WELL_NAME NOT SET");             \
  }                                                                     \
  template <typename strategy_t>                                        \
  BS_API_PLUGIN std::string                                             \
  event_name <strategy_t>::get_event_name () const                      \
  {                                                                     \
    return BOOST_PP_STRINGIZE (event_name);                             \
  }

#define DECL_GET_WELL_EVENT_NAME(event_name)                            \
  template <typename strategy_t>                                        \
  BS_API_PLUGIN std::string                                             \
  event_name <strategy_t>::get_event_name () const                      \
  {                                                                     \
    return BOOST_PP_STRINGIZE (event_name);                             \
  }
  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WELSPECS_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &/*calc_model*/) const
  {
    const sp_top &locked (top);
    // TODO: FIELD
    std::string group_name  = main_params_->get_WELL_GROUP_NAME ("");
    std::string name        = main_params_->get_WELL_NAME ("");

    BS_ASSERT (name.size ());

    sp_well_t well = locked->get_well (group_name, name);
    if (!well)
      {
        well = locked->create_well (group_name, name);
        BS_ERROR (well, "apply_internal");// (group_name) (name);
      }

    const sp_well_t &locked_well (well);

    BS_ASSERT (main_params_->get_I (0) && main_params_->get_J (0)) (main_params_->get_I (0)) (main_params_->get_J (0));

    locked_well->set_coord (main_params_->get_I (), main_params_->get_J ());
    locked_well->set_bhp_depth (main_params_->get_BHP_DEPTH ());
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WELLCON_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &/*calc_model*/) const
  {
    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_GROUP_NAME (""), main_params_->get_WELL_NAME (""));
    if (!well)
      {
        well = locked->create_well (main_params_->get_WELL_GROUP_NAME (""), main_params_->get_WELL_NAME (""));
        if (!well)
          {
            bs_throw_exception (boost::format ("Can't create well (name: %s, group: %s)") % main_params_->get_WELL_NAME ("") % main_params_->get_WELL_GROUP_NAME (""));
          }
      }

    bs_throw_exception ("NOT IMPL YET");

    //const sp_well_t &locked_well (well);
    //locked_well->set_coord (main_params_->get_I (0), main_params_->get_J (0));
    //locked_well->set_state (get_int_d (WELL_STATE));
    //locked_well->set_layer_interval (get_int_d (FIRST_LAYER), get_int_d (LAST_LAYER));
    //locked_well->set_diameter (get_float_d (WELL_DIAMETER));
    //locked_well->set_skin (main_params_->get_SKIN (0));
    //locked_well->set_Kh (get_float_d (KH));
    //locked_well->set_direction (get_int_d (DIRECTION));
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  COMPDAT_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &msh, const sp_calc_model_t &calc_model) const
  {
    const sp_top &locked (top);
    const sp_smesh_iface_t struct_msh (msh, bs_dynamic_cast());
    BS_ASSERT (struct_msh);
    
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    index_t i_coord     = main_params_->get_I (0);
    index_t j_coord     = main_params_->get_J (0);
    index_t first_layer = main_params_->get_FIRST_LAYER (0);
    index_t last_layer  = main_params_->get_LAST_LAYER (0);
    BS_ASSERT (last_layer - first_layer >= 0) (last_layer) (first_layer);

    const sp_well_t &locked_well (well);

    using namespace wells;
    connection_status_type perf_status  = wells::connection_status_cast (main_params_->get_PERFORATION_STATUS (""));
    item_t perf_factor                  = main_params_->get_PERFORATION_FACTOR (0);
    item_t well_diameter                = main_params_->get_WELL_DIAMETER (0);
    item_t kh                           = main_params_->get_KH (0);
    item_t skin                         = main_params_->get_SKIN (0);
    connection_direction_type dir       = wells::connection_direction_cast (main_params_->get_DIRECTION (""));
    index_t seg_number                  = main_params_->get_SEG_NUMBER (0) - 1;
    for (index_t k_coord = first_layer; k_coord != last_layer + 1; ++k_coord)
      {
        BS_ASSERT (i_coord > 0) (i_coord);
        BS_ASSERT (j_coord > 0) (j_coord);
        BS_ASSERT (k_coord > 0) (k_coord);

        index_t n_block = struct_msh->get_element_ijk_to_int (i_coord - 1, j_coord - 1, k_coord - 1);
        BS_ASSERT (n_block >= 0) (i_coord) (j_coord) (k_coord);

        if (n_block < 0)
          {
            BOSOUT (section::schedule, level::warning) << boost::format ("COMPDAT: Connection for well %s in inactive cell [%d, %d, %d]") % well->name () % i_coord % j_coord % k_coord << bs_end;
            continue;
          }
        else
          {
            BOSOUT (section::schedule, level::debug) << boost::format ("COMPDAT: Connection for well %s in cell [%d, %d, %d, %d]") % well->name () % i_coord % j_coord % k_coord % n_block << bs_end;
          }

        sp_connection_t con = well->get_connection_map (n_block);
        if (!con)
          {
            con = well->add_primary_connection (i_coord - 1, j_coord - 1, k_coord - 1, n_block);
          }
        else
          {
            BOSOUT (section::schedule, level::debug) << boost::format ("COMPDAT: Connection for well %s in cell [%d, %d, %d, %d] will be set to [%d, %d, %d, %d]") 
              % well->name () % (con->i_coord () + 1) % (con->j_coord () + 1) % (con->k_coord () + 1) % con->n_block ()
              % i_coord % j_coord % k_coord % n_block 
              << bs_end;
          }

        const sp_connection_t &locked_connection (con);

        locked_connection->set_status (perf_status);
        locked_connection->set_factor (perf_factor);
        locked_connection->set_diameter (well_diameter);
        locked_connection->set_Kh (kh);
        locked_connection->set_skin (skin);
        locked_connection->set_direction (dir);
        locked_connection->set_connection_depth (msh);
        locked_connection->set_mult (1.0);
        locked_connection->set_seg_number (seg_number);
      }

    locked_well->check_shut ();
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WCONPROD_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &calc_model) const
  {
    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    const sp_well_t &locked_well (well);
    locked_well->set_state (well_state_cast (main_params_->get_WELL_STATE ("")), calc_model);

    sp_well_controller_t controller = well->get_controller ();
    if (!controller)
      controller = locked->create_well_controller (well);

    const sp_well_controller_t &locked_controller (controller);
    BS_ASSERT (locked_controller);

    using namespace wells;

    locked_controller->set_bhp (main_params_->get_BHP (calc_model->internal_constants.default_production_bhp_limit));

    // TODO: BUG:
    locked_well->set_bhp (locked_controller->bhp ());

    wells::rate_control_type control_type = wells::rate_control_cast (main_params_->get_WELL_CONTROL_MODE (""));
    locked_controller->set_main_control (well, control_type, true);

    // TODO: bad design, but i don't know better way
    locked_controller->clear_rate ();
    locked_controller->set_rate (oil_rate_value, main_params_->get_OIL_RATE (0));
    locked_controller->set_rate (water_rate_value, main_params_->get_WATER_RATE (0));
    locked_controller->set_rate (gas_rate_value, main_params_->get_GAS_RATE (0));
    locked_controller->set_rate (liquid_rate_value, main_params_->get_OUTER_RATE (calc_model->internal_constants.default_liquid_rate_limit));
    locked_controller->set_rate (liquid_inner_rate_value, main_params_->get_INNER_RATE (0));
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WCONHIST_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &calc_model) const
  {
    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    const sp_well_t &locked_well (well);
    locked_well->set_state (well_state_cast (main_params_->get_WELL_STATE ("")), calc_model);

    sp_well_controller_t controller = well->get_controller ();
    if (!controller)
      controller = locked->create_well_controller (well);

    const sp_well_controller_t &locked_controller (controller);
    BS_ASSERT (locked_controller);

    using namespace wells;
    locked_controller->set_bhp (calc_model->internal_constants.default_production_bhp_limit);
    locked_controller->set_bhp_history (main_params_->get_BHP (calc_model->internal_constants.default_production_bhp_limit));

    //BOSOUT (section::schedule, level::debug) << "default_production_bhp_limit: " << calc_model->internal_constants.default_production_bhp_limit << bs_end;

    // TODO: BUG:
    locked_well->set_bhp (locked_controller->bhp ());

    wells::rate_control_type control_type = wells::rate_control_cast (main_params_->get_WELL_CONTROL_MODE (""));
    locked_controller->set_main_control (well, control_type, true);

    // TODO: bad design, but i don't know better way
    locked_controller->clear_rate ();
    locked_controller->set_rate (liquid_rate_value, calc_model->internal_constants.default_liquid_rate_limit);
    locked_controller->set_rate (oil_rate_value, main_params_->get_OIL_RATE (0));
    locked_controller->set_rate (water_rate_value, main_params_->get_WATER_RATE (0));
    locked_controller->set_rate (gas_rate_value, main_params_->get_GAS_RATE (0));
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WCONINJE_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &calc_model) const
  {
    using namespace wells;

    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    const sp_well_t &locked_well (well);
    locked_well->set_state (well_state_cast (main_params_->get_WELL_STATE ("")), calc_model);

    sp_well_controller_t controller = well->get_controller ();
    if (!controller)
      controller = locked->create_well_controller (well);

    const sp_well_controller_t &locked_controller (controller);
    BS_ASSERT (locked_controller);

    locked_controller->set_injection_type (wells::injection_type_cast (main_params_->get_INJECTOR_TYPE ("NONE")));
    locked_controller->set_bhp (main_params_->get_BHP (calc_model->internal_constants.default_injection_bhp_limit));

    // TODO: BUG:
    locked_well->set_bhp (locked_controller->bhp ());

    wells::rate_control_type control_type = wells::rate_control_cast (main_params_->get_WELL_CONTROL_MODE (""));
    locked_controller->set_main_control (well, control_type, false);

    // TODO: bad design, but i don't know better way
    locked_controller->clear_rate ();
    locked_controller->set_rate (rate_value, main_params_->get_OUTER_RATE (calc_model->internal_constants.default_liquid_rate_limit));
    locked_controller->set_rate (liquid_inner_rate_value, main_params_->get_INNER_RATE (0));
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WECON_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &/*calc_model*/) const
  {
    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    wells::limit_operation_type type = wells::limit_operation_cast (main_params_->get_OPERATION (""));
    const sp_well_limit_operation_t &locked_limit (locked->create_well_limit_operation (well, type));
    BS_ASSERT (locked_limit) (main_params_->get_WELL_NAME ("")) (main_params_->get_OPERATION (""));

    locked_limit->set_min_oil_rate (main_params_->get_MIN_OIL_RATE (0));
    locked_limit->set_max_water_cut (main_params_->get_WATER_CUT (0));
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WECONINJ_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &/*calc_model*/) const
  {
    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    if (!well->get_limit_operation ())
      locked->create_well_limit_operation (well, wells::operation_none);

    const sp_well_limit_operation_t &locked_limit (well->get_limit_operation ());
    locked_limit->set_min_water_cut (main_params_->get_MIN_RATE (0));
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WEFAC_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &/*calc_model*/) const
  {
    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    const sp_well_t &locked_well (well);
    locked_well->set_exploitation_factor (main_params_->get_OPERATION_FACTOR (-1.0));
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WELTARG_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &calc_model) const
  {
    const sp_top &locked (top);
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    BS_ASSERT (well->get_controller ());
    const sp_well_controller_t &locked_controller (well->get_controller ());

    using namespace wells;
    rate_value_type type = wells::rate_value_cast (main_params_->get_WELL_CONTROL (""));
    if (type == bhp_value)
      {
        typename strategy_t::item_t bhp = locked_controller->is_production ()
          ? main_params_->get_VALUE (calc_model->internal_constants.default_production_bhp_limit)
          : main_params_->get_VALUE (calc_model->internal_constants.default_injection_bhp_limit)
          ;

        if (!well->is_rate ())
          well->set_bhp (bhp);

        locked_controller->set_bhp (bhp);
      }
    else
      {
        locked_controller->set_rate (type, main_params_->get_VALUE (0));
      }
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  WPIMULT_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &msh, const sp_calc_model_t &/*calc_model*/) const
  {
    const sp_top &locked (top);
    const sp_smesh_iface_t struct_msh (msh, bs_dynamic_cast());
    
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    typedef typename strategy_t::index_t index_t;

    double perm_mult = main_params_->get_PERM_FACTOR (0);
    index_t i_cell   = main_params_->get_I (0);
    index_t j_cell   = main_params_->get_J (0);
    index_t k_cell   = main_params_->get_K (0);
    index_t z1_cell  = main_params_->get_Z1 (0);
    index_t z2_cell  = main_params_->get_Z2 (0);

    if (i_cell && j_cell && k_cell && z1_cell && z2_cell)
      {
        bs_throw_exception ("Should be specified only [i, j, k] or [z1, z2] values");
      }
    if (z1_cell || z2_cell)
      {
        bs_throw_exception ("[z1, z2] values not supported now");
      }

    if (/*(!z1_cell && !z2_cell) || */!i_cell || !j_cell || !k_cell)
      {
        typename well_t::connection_iterator_t it = well->connections_begin (),
                 e = well->connections_end ();
        for (; it != e; ++it)
          {
            it->mul_perm_mult (perm_mult);
          }
      }
    else if (i_cell && j_cell && k_cell)
      {
        index_t n_block = struct_msh->get_element_ijk_to_int (i_cell - 1, j_cell - 1, k_cell - 1);
        const typename well_t::sp_connection_t &c = well->get_connection_map (n_block);
        if (c)
          {
            c->mul_perm_mult (perm_mult);
          }
      }
    //else 
    //  {
    //    if (!z1_cell)
    //      {
    //        z1_cell = 1;
    //      }

    //    if (!z2_cell)
    //      {
    //        z2_cell = struct_msh->get_dimens()[2];
    //      }

    //    index_t x = well->i_coord_;
    //    index_t y = well->j_coord_;
    //    for (index_t z = z1_cell; z < z2_cell; ++z)
    //      {
    //        index_t n_block = struct_msh->get_element_ijk_to_int (x, y, z - 1);
    //        const typename well_t::sp_connection_t &c = well->get_connection (n_block);
    //        if (c)
    //          {
    //            c->mul_perm_mult (perm_mult);
    //          }
    //      }
    //  }
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  COMPENSATION_event<strategy_t>::apply_internal (const sp_top &/*top*/, const sp_mesh_iface_t &/*msh*/, const sp_calc_model_t &/*calc_model*/) const
  {
    bs_throw_exception ("NOT IMPL YET");
  }

  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  FRACTURE_event<strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &msh, const sp_calc_model_t &calc_model) const
  {
    const sp_top &locked (top);
    const sp_smesh_iface_t struct_msh (msh, bs_dynamic_cast());
    
    sp_well_t well = locked->get_well (main_params_->get_WELL_NAME (""));
    BS_ASSERT (well) (main_params_->get_WELL_NAME (""));
    if (!well)
      return ;

    typedef typename strategy_t::index_t index_t;
    typedef typename strategy_t::item_t  item_t;

    index_t i_cell = main_params_->get_I (0);
    index_t j_cell = main_params_->get_J (0);
    index_t kw1_cell = main_params_->get_KW1 (0);
    index_t kw2_cell = main_params_->get_KW2 (0);

    if (!i_cell || !j_cell || !kw1_cell || !kw2_cell)
      {
        bs_throw_exception ("Default values not allowed");
      }

    item_t half_length  = main_params_->get_HALF_LENGHT (0);
    item_t theta        = main_params_->get_TETTA (0);
    item_t skin         = main_params_->get_SKIN (0);

    for (index_t k = kw1_cell; k < kw2_cell; ++k)
      {
        index_t n_block = struct_msh->get_element_ijk_to_int (i_cell - 1, j_cell - 1, k - 1);
        const typename well_t::sp_connection_t &c = well->get_connection_map (n_block);
        if (c)
          {
            c->set_half_length (half_length);
            c->set_theta (theta);
            c->set_skin (skin);
          }

      }

    well->set_state (well_state_cast (main_params_->get_WELL_STATE ("")), calc_model);
  }

#define SPEC_APPLY(name)																																																		          \
	template void name<base_strategy_fi>::apply_internal (const sp_top &top, const sp_mesh_iface_t &msh, const sp_calc_model_t &calc_model) const;\
	template void name<base_strategy_di>::apply_internal (const sp_top &top, const sp_mesh_iface_t &msh, const sp_calc_model_t &calc_model) const;\
	template void name<base_strategy_mixi>::apply_internal (const sp_top &top, const sp_mesh_iface_t &msh, const sp_calc_model_t &calc_model) const;

  SPEC_APPLY (WELSPECS_event);
  SPEC_APPLY (WELLCON_event);
  SPEC_APPLY (COMPDAT_event);
  SPEC_APPLY (WCONPROD_event);
  SPEC_APPLY (WCONHIST_event);
  SPEC_APPLY (WCONINJE_event);
  SPEC_APPLY (WECON_event);
  SPEC_APPLY (WECONINJ_event);
  SPEC_APPLY (WEFAC_event);
  SPEC_APPLY (WELTARG_event);
  SPEC_APPLY (WPIMULT_event);
  SPEC_APPLY (COMPENSATION_event);
  SPEC_APPLY (FRACTURE_event);
//SPEC_APPLY (PERMFRAC_event);

  DECL_GET_WELL_EVENT_NAMES (WELSPECS_event);
  DECL_GET_WELL_EVENT_NAMES (WELLCON_event);
  DECL_GET_WELL_EVENT_NAMES (COMPDAT_event);
  DECL_GET_WELL_EVENT_NAMES (WCONPROD_event);
  DECL_GET_WELL_EVENT_NAMES (WCONHIST_event);
  DECL_GET_WELL_EVENT_NAMES (WCONINJE_event);
  DECL_GET_WELL_EVENT_NAMES (WECON_event);
  DECL_GET_WELL_EVENT_NAMES (WECONINJ_event);
  DECL_GET_WELL_EVENT_NAMES (WEFAC_event);
  DECL_GET_WELL_EVENT_NAMES (WELTARG_event);
  DECL_GET_WELL_EVENT_NAMES (WPIMULT_event);
  DECL_GET_WELL_EVENT_NAMES (FRACTURE_event);

  DECL_GET_WELL_EVENT_NAME (COMPENSATION_event);
  DECL_GET_WELL_EVENT_NAME (PERMFRAC_event);

  template class WELSPECS_event <base_strategy_fi>;
  template class WELSPECS_event <base_strategy_di>;
  template class WELSPECS_event <base_strategy_mixi>;
  template class WELLCON_event <base_strategy_fi>;
  template class WELLCON_event <base_strategy_di>;
  template class WELLCON_event <base_strategy_mixi>;
  template class COMPDAT_event <base_strategy_fi>;
  template class COMPDAT_event <base_strategy_di>;
  template class COMPDAT_event <base_strategy_mixi>;
  template class WCONPROD_event <base_strategy_fi>;
  template class WCONPROD_event <base_strategy_di>;
  template class WCONPROD_event <base_strategy_mixi>;
  template class WCONHIST_event <base_strategy_fi>;
  template class WCONHIST_event <base_strategy_di>;
  template class WCONHIST_event <base_strategy_mixi>;
  template class WCONINJE_event <base_strategy_fi>;
  template class WCONINJE_event <base_strategy_di>;
  template class WCONINJE_event <base_strategy_mixi>;
  template class WECON_event <base_strategy_fi>;
  template class WECON_event <base_strategy_di>;
  template class WECON_event <base_strategy_mixi>;
  template class WECONINJ_event <base_strategy_fi>;
  template class WECONINJ_event <base_strategy_di>;
  template class WECONINJ_event <base_strategy_mixi>;
  template class WEFAC_event <base_strategy_fi>;
  template class WEFAC_event <base_strategy_di>;
  template class WEFAC_event <base_strategy_mixi>;
  template class WELTARG_event <base_strategy_fi>;
  template class WELTARG_event <base_strategy_di>;
  template class WELTARG_event <base_strategy_mixi>;
  template class WPIMULT_event <base_strategy_fi>;
  template class WPIMULT_event <base_strategy_di>;
  template class WPIMULT_event <base_strategy_mixi>;
  template class COMPENSATION_event <base_strategy_fi>;
  template class COMPENSATION_event <base_strategy_di>;
  template class COMPENSATION_event <base_strategy_mixi>;
  template class FRACTURE_event <base_strategy_fi>;
  template class FRACTURE_event <base_strategy_di>;
  template class FRACTURE_event <base_strategy_mixi>;
  template class PERMFRAC_event <base_strategy_fi>;
  template class PERMFRAC_event <base_strategy_di>;
  template class PERMFRAC_event <base_strategy_mixi>;

} // namespace blue_sky

