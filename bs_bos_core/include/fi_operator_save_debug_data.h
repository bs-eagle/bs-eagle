/**
 * \file fi_operator_save_debug_data.h
 * \brief save debug data on each fi_operator iteration
 * \author Sergey Miryanov
 * \date 03.02.2009
 * */
#ifndef BS_FI_OPERATOR_SAVE_DEBUG_DATA_H_
#define BS_FI_OPERATOR_SAVE_DEBUG_DATA_H_

#include "save_connection_data.h"
#include "string_formater.h"
#include "calc_well.h"
#include "well_iterator.h"
#include "rr_rw_wr_saver.h"
#include "well_type_helper.h"
#include "member_accessor.h"

#define SAVE_BOOST_ARRAY(name, size_j) \
	tools::save_seq_vector (tools::string_formater (std::string (BOOST_PP_STRINGIZE(name)) + std::string (".bs.%d.txt"), iter_counter).str) \
		.save_via_fn (debug::BOOST_PP_CAT (save_, name) <data_array_t, item_t> (calc_model_->data, size_j));

#define SAVE_ITEM(name) \
	tools::save_seq_vector (tools::string_formater (std::string (BOOST_PP_STRINGIZE(name)) + std::string (".bs.%d.txt"), iter_counter).str) \
		.save_via_fn (debug::BOOST_PP_CAT (save_, name) <data_array_t, item_t> (calc_model_->data));

namespace blue_sky {

  namespace debug {

    template <typename well_t>
    struct ww_diveder
    {
      typedef typename well_t::item_t item_t;
      item_t operator () (const well_t *well) const
      {
        return 1.0 / well->get_ww_value ();
      }
    };
    template <typename well_t>
    struct ww_accessor
    {
      typedef typename well_t::item_t item_t;
      item_t operator () (const well_t *well) const
      {
        return well->get_ww_value ();
      }
    };
    template <typename well_t>
    struct bw_accessor
    {
      typedef typename well_t::item_t item_t;
      item_t operator () (const well_t *well) const
      {
        return well->get_bw_value ();
      }
    };
    
    template <typename connection_t>
    struct rr
    {
      typedef typename connection_t::item_t item_t;
      item_t operator () (const connection_t *c, size_t i) const
      {
        return c->get_rr_value () [i];
      };
    };
    template <typename connection_t>
    struct rw
    {
      typedef typename connection_t::item_t item_t;
      item_t operator () (const connection_t *c, size_t i) const
      {
        return c->get_rw_value () [i];
      };
    };
    template <typename connection_t>
    struct wr
    {
      typedef typename connection_t::item_t item_t;
      item_t operator () (const connection_t *c, size_t i) const
      {
        return c->get_wr_value () [i];
      };
    };
    template <typename connection_t>
    struct rate 
    { 
      typedef typename connection_t::item_t item_t;
      item_t operator () (const connection_t *c, size_t i) const 
      { 
        return c->get_rate_value () [i]; 
      }
    };
    
    template <typename connection_t>
    struct cur_bhp
    {
      typedef typename connection_t::item_t item_t;
      item_t operator () (const connection_t *c) const
      {
        return c->get_cur_bhp ();
      };
    };
    template <typename connection_t>
    struct density
    {
      typedef typename connection_t::item_t item_t;
      item_t operator () (const connection_t *c) const
      {
        return c->density;
      };
    };
    template <typename connection_t>
    struct bulkp
    {
      typedef typename connection_t::item_t item_t;
      item_t operator () (const connection_t *c) const
      {
        return c->bulkp;
      };
    };

#define MEMBER_SAVER(name)                                \
    template <typename data_array_t, typename item_t>     \
    struct BOOST_PP_CAT (save_, name)                     \
    {                                                     \
      typedef item_t value_type;                          \
                                                          \
      BOOST_PP_CAT (save_, name) (const data_array_t &data_, size_t size_j_)  \
      : data_ (data_)                                     \
      , size_j_ (size_j_)                                 \
      {                                                   \
      }                                                   \
                                                          \
      size_t                                              \
      size_i () const                                     \
      {                                                   \
        return data_.size ();                             \
      }                                                   \
                                                          \
      size_t                                              \
      size_j () const                                     \
      {                                                   \
        return size_j_;                                   \
      }                                                   \
                                                          \
      item_t                                              \
      get (size_t i, size_t j) const                      \
      {                                                   \
        return data_[i].name[j];                          \
      }                                                   \
                                                          \
      const data_array_t &data_;                          \
      size_t size_j_;                                     \
    };

#define MEMBER_SAVER_2(name)                              \
    template <typename data_array_t, typename item_t>     \
    struct BOOST_PP_CAT (save_, name)                     \
    {                                                     \
      typedef item_t value_type;                          \
                                                          \
      BOOST_PP_CAT (save_, name) (const data_array_t &data_)  \
      : data_ (data_)                                     \
      {                                                   \
      }                                                   \
                                                          \
      size_t                                              \
      size_i () const                                     \
      {                                                   \
        return data_.size ();                             \
      }                                                   \
                                                          \
      size_t                                              \
      size_j () const                                     \
      {                                                   \
        return 1;                                         \
      }                                                   \
                                                          \
      item_t                                              \
      get (size_t i, size_t j) const                      \
      {                                                   \
        return data_[i].name;                             \
      }                                                   \
                                                          \
      const data_array_t &data_;                          \
    };

    MEMBER_SAVER (mobility);
    MEMBER_SAVER (p_deriv_mobility);
    MEMBER_SAVER (s_deriv_mobility);
    MEMBER_SAVER (density);
    MEMBER_SAVER (p_deriv_density);
    MEMBER_SAVER (relative_perm);
    MEMBER_SAVER (s_deriv_relative_perm);
    MEMBER_SAVER (invers_visc_fvf);
    MEMBER_SAVER (p_deriv_invers_visc_fvf);
    MEMBER_SAVER (cap_pressure);
    MEMBER_SAVER (s_deriv_cap_pressure);

    MEMBER_SAVER_2 (gor_deriv_invers_fvf);
    MEMBER_SAVER_2 (gor_deriv_invers_viscosity);
    MEMBER_SAVER_2 (gor_deriv_invers_visc_fvf);

  }
  /**
   * \brief debug function for output data
   * \param dt current dt step
   * */
  template <typename strategy_t, bool is_w, bool is_g, bool is_o>
  void
  fi_operator_impl <strategy_t, is_w, is_g, is_o>::debug_save_data (item_t dt)
  {
    static int iter_counter = 0;
    iter_counter++;
    //int zxczxc = iter_counter;

    //tools::save_seq_vector (tools::string_formater ("rhs.bs.%d.txt", iter_counter).str).save (jmatrix_->get_rhs ());
    //tools::save_seq_vector (tools::string_formater ("rhs_flux.bs.%d.txt", iter_counter).str).save (jmatrix_->get_rhs_flux ());
    //tools::save_seq_vector (tools::string_formater ("sol.bs.%d.txt", iter_counter).str).save (jmatrix_->get_solution ());
    //tools::save_seq_vector (tools::string_formater ("regular_diag.bs.%d.txt", iter_counter).str).save (jmatrix_->get_regular_acc_diag ());

    //jmatrix_->get_regular_matrix ()->ascii_write_in_csr_format (tools::string_formater ("regular_csr_matrix.bs.%d.txt", iter_counter).str);
    //jmatrix_->get_irregular_matrix ()->ascii_write_in_csr_format (tools::string_formater ("irregular_csr_matrix.bs.%d.txt", iter_counter).str);

    //SAVE_BOOST_ARRAY (cap_pressure, calc_model_->n_phases - 1);
    //SAVE_BOOST_ARRAY (s_deriv_cap_pressure, calc_model_->n_phases - 1);
    //SAVE_BOOST_ARRAY (relative_perm, calc_model_->n_phases);
    //SAVE_BOOST_ARRAY (s_deriv_relative_perm, (2 * (calc_model_->n_phases - 1)));
    //////tools::save_seq_vector (tools::string_formater ("gor.bs.%d.txt", iter_counter).str).save (calc_model_->gas_oil_ratio);
    //////SAVE_ITEM (p_deriv_gas_oil_ratio);

    //SAVE_BOOST_ARRAY (invers_fvf);
    //SAVE_BOOST_ARRAY (p_deriv_invers_fvf);
    //SAVE_ITEM (gor_deriv_invers_fvf);

    //SAVE_BOOST_ARRAY (invers_viscosity);
    //SAVE_BOOST_ARRAY (p_deriv_invers_viscosity);
    //SAVE_ITEM (gor_deriv_invers_viscosity);

    //SAVE_BOOST_ARRAY (invers_visc_fvf, calc_model_->n_phases);
    //SAVE_BOOST_ARRAY (p_deriv_invers_visc_fvf, calc_model_->n_phases);
    //SAVE_ITEM (gor_deriv_invers_visc_fvf);

    ////SAVE_ITEM (porosity);
    ////SAVE_ITEM (p_deriv_porosity);

    //SAVE_BOOST_ARRAY (density, calc_model_->n_phases);
    //SAVE_BOOST_ARRAY (p_deriv_density, calc_model_->n_phases);

    //SAVE_BOOST_ARRAY (mobility, calc_model_->n_phases);
    //SAVE_BOOST_ARRAY (p_deriv_mobility, calc_model_->n_phases);
    //SAVE_BOOST_ARRAY (s_deriv_mobility, (2 * (calc_model_->n_phases - 1)));

    //tools::save_seq_vector (tools::string_formater ("sat_3p.bs.%d.txt", iter_counter).str).save (calc_model_->saturation_3p);
    //tools::save_seq_vector (tools::string_formater ("pressure.bs.%d.txt", iter_counter).str).save (calc_model_->pressure);
    //tools::save_seq_vector (tools::string_formater ("gas_oil_ratio.bs.%d.txt", iter_counter).str).save (calc_model_->gas_oil_ratio);

    //tools::save_seq_vector ("perm.bs.txt").save (tools::float_saver (calc_model_->rock_grid_prop->permeability));
    //tools::save_seq_vector ("poro.bs.txt").save (tools::float_saver (calc_model_->rock_grid_prop->porosity_p_ref));
    //tools::save_seq_vector ("ngt.bs.txt").save (tools::float_saver (calc_model_->rock_grid_prop->net_to_gros));
    //tools::save_seq_vector ("volume.bs.txt").save (tools::float_saver (calc_model_->rock_grid_prop->volume));

    //tools::save_seq_vector (tools::string_formater ("s_rhs.bs.%d.txt", iter_counter).str).save (jmatrix_->get_sec_rhs ());
    //tools::save_seq_vector (tools::string_formater ("ss_diag.bs.%d.txt", iter_counter).str).save (jmatrix_->get_ss_diagonal ());
    //tools::save_seq_vector (tools::string_formater ("sp_diag.bs.%d.txt", iter_counter).str).save (jmatrix_->get_sp_diagonal ());

    typedef typename calc_model_t::reservoir_t::facility_manager_t::well_const_iterator_t well_iterator_t;
    well_iterator_t wb = reservoir_->get_facility_list ()->wells_begin ();
    well_iterator_t we = reservoir_->get_facility_list ()->wells_end ();

    typedef wells::type_helper <strategy_t>          type_helper_t;
    typedef typename type_helper_t::item_rr_block_t item_rr_block_t;
    typedef typename type_helper_t::item_rw_block_t item_rw_block_t;
    typedef typename type_helper_t::item_wr_block_t item_wr_block_t;

    using namespace debug;
    typedef rr          <wells::connection <strategy_t> > rr_t;
    typedef rw          <wells::connection <strategy_t> > rw_t;
    typedef wr          <wells::connection <strategy_t> > wr_t;
    typedef rate        <wells::connection <strategy_t> > rate_t;
    typedef cur_bhp     <wells::connection <strategy_t> > cur_bhp_t;
    typedef density     <wells::connection <strategy_t> > density_t;
    typedef bulkp       <wells::connection <strategy_t> > bulkp_t;
    typedef ww_diveder  <well <strategy_t> >              ww_diveder_t;
    typedef bw_accessor <well <strategy_t> >              bw_accessor_t;

    //tools::save_seq_vector (tools::string_formater ("rr.bs.%d.txt", iter_counter).str)
    //.save (tools::rr_rw_wr_saver <strategy_t, item_rr_block_t, rr_t> (wb, we, reservoir_->get_connections_count ()));
    //tools::save_seq_vector (tools::string_formater ("rw.bs.%d.txt", iter_counter).str)
    //.save (tools::rr_rw_wr_saver <strategy_t, item_rw_block_t, rw_t> (wb, we, reservoir_->get_connections_count ()));
    //tools::save_seq_vector (tools::string_formater ("wr.bs.%d.txt", iter_counter).str)
    //.save (tools::rr_rw_wr_saver <strategy_t, item_wr_block_t, wr_t> (wb, we, reservoir_->get_connections_count ()));
    //tools::save_seq_vector (tools::string_formater ("rate.bs.%d.txt", iter_counter).str)
    //.save (tools::rr_rw_wr_saver <strategy_t, item_wr_block_t, rate_t> (wb, we, reservoir_->get_connections_count ()));

    //tools::save_seq_vector (tools::string_formater ("con_cur_bhp.bs.%d.txt", iter_counter).str)
    //.save (tools::connection_member_saver <strategy_t, cur_bhp_t> (wb, we, reservoir_->get_connections_count ()));
    //tools::save_seq_vector (tools::string_formater ("con_density.bs.%d.txt", iter_counter).str)
    //.save (tools::connection_member_saver <strategy_t, density_t> (wb, we, reservoir_->get_connections_count ()));
    //tools::save_seq_vector (tools::string_formater ("con_bulkp.bs.%d.txt", iter_counter).str)
    //.save (tools::connection_member_saver <strategy_t, bulkp_t> (wb, we, reservoir_->get_connections_count ()));

    //tools::save_seq_vector (tools::string_formater ("ww.bs.%d.txt", iter_counter).str)
    //.save (tools::well_member_saver <strategy_t, ww_diveder_t> (wb, we));
    //tools::save_seq_vector (tools::string_formater ("bw.bs.%d.txt", iter_counter).str)
    //.save (tools::well_member_saver <strategy_t, bw_accessor_t> (wb, we));

    //static save_connection_data <strategy_t> connection_data;
    //connection_data.save ("con_data.bs.%d.txt", calc_model_, dt, reservoir_->get_facility_list ()->wells_begin (), reservoir_->get_facility_list ()->wells_end (), iter_counter, 0);

    BOSOUT (section::iters, level::debug) << "iter_counter: " << iter_counter << bs_end;
  }

} // namespace blue_sky


#endif  // #ifndef BS_FI_OPERATOR_SAVE_DEBUG_DATA_H_

