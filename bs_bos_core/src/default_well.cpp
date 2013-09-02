/**
 *       \file  default_well.cpp
 *      \brief  Implementation of default_well
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  20.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "default_well.h"
#include "default_rr_eliminator.h"
#include "default_connection.h"
#include "calc_model.h"

// TODO: 
#include "default_well_calc_rate.h"
#include "default_connection_iterator.h"

#include "add_connection_to_list.h"
#include "well_tools.h"

namespace blue_sky {
namespace wells {

  /**
   * \brief  'default' ctor for default_well
   * \param  param Additional params for ctor
   * */
  default_well::default_well (bs_type_ctor_param param /* = NULL */)
  : base_t (param)
  , ww_value (0)
  , bw_value (0)
  , open_connections_count_ (0)
  {
  }

  /**
   * \brief  ctor for default_well
   * \param  well_name Name of well
   * \todo   Obsolete 
   * */
  default_well::default_well (const std::string &well_name)
  : base_t (well_name)
  , ww_value (0)
  , bw_value (0)
  , open_connections_count_ (0)
  {

  }

  /**
   * \brief  copy-ctor for default_well
   * \param  w Instance of default_well to be copied
   * */
  default_well::default_well (const default_well &w)
        : bs_refcounter (), well ()
  {
    *this = w;
  }

  default_well::sp_connection_t
    default_well::add_completion (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block, completion_coords_t &completion_coords_data)
  {
    sp_connection_t con = detail::add_connection (i_coord, j_coord, k_coord, n_block, primary_connection_list_, this);
    con->set_completion_coords (completion_coords_data);
    return con;
  }

  default_well::sp_connection_t
  default_well::add_primary_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block)
  {
    return detail::add_connection (i_coord, j_coord, k_coord, n_block, primary_connection_list_, this);
  }
  default_well::sp_connection_t
  default_well::add_secondary_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block)
  {
    return detail::add_connection (i_coord, j_coord, k_coord, n_block, secondary_connection_list_, this);
  }

  default_well::sp_connection_t
  default_well::get_connection_map (index_t n_block) const
  {
    for (size_t i = 0, cnt = primary_connection_list_.size (); i < cnt; ++i)
      {
        if (primary_connection_list_[i]->n_block () == n_block)
          return primary_connection_list_[i];
      }

    for (size_t i = 0, cnt = secondary_connection_list_.size (); i < cnt; ++i)
      {
        if (secondary_connection_list_[i]->n_block () == n_block)
          return secondary_connection_list_[i];
      }

    return 0;
//    typename connection_map_t::const_iterator it = connection_map_.find (n_block);
//    if (it == connection_map_.end ())
//      return 0;
//
//    if (it->second >= (int)connection_list_.size ())
//      {
//#ifdef _DEBUG
//        typename connection_map_t::const_iterator i = connection_map_.begin (), e = connection_map_.end ();
//        for (; i != e; ++i)
//          {
//            BOSOUT (section::wells, level::debug) << boost::format ("[%d: %d]") % i->first % i->second << bs_end;
//          }
//#endif
//        bs_throw_exception (boost::format ("Connection map has an invalid value ([%s]: %d - %d, %d, %d)") % base_t::name_ % n_block % it->second % connection_list_.size () % connection_map_.size ());
//      }
//
//    return connection_list_[it->second];
  }

  bool
  default_well::is_primary_connection (index_t n_block) const
  {
    for (size_t i = 0, cnt = primary_connection_list_.size (); i < cnt; ++i)
      {
        if (primary_connection_list_[i]->n_block () == n_block)
          return true;
      }

    return false;
  }

  bool
  default_well::is_primary_connection (const connection_iterator &it) const
  {
    return it.position () < static_cast <index_t> (primary_connection_list_.size ());
  }

  connection_iterator
  default_well::connections_begin () const
  {
    return base_t::connection_iterator_t (
      new default_connection_iterator <default_well, default_connection_t> (this, begin_iterator_tag));
  }

  connection_iterator
  default_well::connections_end () const
  {
    return base_t::connection_iterator_t (
      new default_connection_iterator <default_well, default_connection_t> (this, end_iterator_tag));
  }

  bool
  default_well::is_no_connections () const
  {
    return primary_connection_list_.empty () && 
      secondary_connection_list_.empty ();
  }

  bool
  default_well::is_no_primary_connections () const
  {
    return primary_connection_list_.empty ();
  }

  default_well::sp_connection_t
  default_well::get_first_connection () const
  {
    if (primary_connection_list_.empty ())
      {
        bs_throw_exception (boost::format ("No primary connections (well: %s)") % base_t::name ());
      }

    return primary_connection_list_[0];
  }


  void
  default_well::restore_solution (double /*dt*/, const spv_double &p_sol, 
                                               const spv_double & /*s_sol*/, index_t block_size)
  {
    BS_ASSERT (!is_no_connections ()) (base_t::name ());

    typedef shared_vector <item_t> item_wr_block_t;
    typedef shared_vector <item_t> item_xr_block_t;

    // TODO: may be we can not check is_work flag
    if (!base_t::well_state_.is_work)
      return ;

    if (!base_t::well_controller_->is_rate ())
      return ;

    item_t xw = 0;
    item_t ww_bw = (1.0 / ww_value) * bw_value;

    typedef default_connection_iterator_impl <default_well, default_connection> iterator_t;
    iterator_t it (this, begin_iterator_tag), e (this, end_iterator_tag);
    for (; !base_t::is_shut () && it != e; ++it)
      {
        const sp_connection_t &c  = *it;
        if (!c->is_shut ())
          {
            const item_wr_block_t &wr = c->get_wr_value ();
            // FIXME: remove shared_array and vectors
            const item_xr_block_t &xr = shared_array (const_cast <item_t *> (&(*p_sol)[c->n_block () * block_size]), block_size);

            item_t wr_xr = 0;
            for (index_t j = 0; j < block_size; ++j)
              {
                wr_xr += wr[j] * xr[j];
              }

            xw += ww_bw - wr_xr;
          }
      }

    BOSOUT (section::wells, level::low) << "[" << base_t::name_ << boost::format ("] Restore solution: change bhp from %f to %f (bw: %f, ww: %f)") % base_t::bhp_ % (base_t::bhp_ + xw) % bw_value % ww_value << bs_end;
    base_t::bhp_ += xw;
  }

  /**
   * \brief  Calculates Jacobian value and stores it in array
   * \param  well
   * \param  array Array of Jacobian values
   * \param  rw_index
   * \param  wr_index
   * \param  dt
   * \param  block_size
   * */
  template <typename well_t, typename rhs_item_t, typename index_t, typename connection_t>
  static void
  eliminate (well_t *well, 
    rhs_item_t *res, 
    const smart_ptr <connection_t> &rw_con, 
    const smart_ptr <connection_t> &wr_con, 
    double dt, 
    index_t block_size) 
  {
    BS_ASSERT (!well->is_shut ()) (well->name ());
    BS_ASSERT (res);

    if (block_size > 3 && block_size < 1)
      bs_throw_exception ("block_size invalid");

    typedef typename well_t::item_t item_t;

    const boost::array <item_t, FI_PHASE_TOT> &rw                 = rw_con->rw_value;
    const boost::array <item_t, FI_PHASE_TOT> &wr                 = wr_con->wr_value;
    const boost::array <item_t, connection_t::rr_value_count> &rr = wr_con->rr_value;

    BS_ASSERT (rw.size () >= (size_t)block_size) (rw.size ()) (block_size);
    BS_ASSERT (wr.size () >= (size_t)block_size) (wr.size ()) (block_size);
    BS_ASSERT (rr.size () >= (size_t)block_size * block_size) (rr.size ()) (block_size);

    if (rw_con->n_block () == wr_con->n_block ())
      {
        if (well->get_well_controller ()->is_rate ())
          {
            if (block_size == 3)
              default_rr_eliminator <3>::process_diag_rate (res, rr, rw, wr, dt);
            else if (block_size == 2)
              default_rr_eliminator <2>::process_diag_rate (res, rr, rw, wr, dt);
            else if (block_size == 1)
              default_rr_eliminator <1>::process_diag_rate (res, rr, rw, wr, dt);
          }
        else
          {
            if (block_size == 3)
              default_rr_eliminator <3>::process_diag_bhp (res, rr, dt);
            else if (block_size == 2)
              default_rr_eliminator <2>::process_diag_bhp (res, rr, dt);
            else if (block_size == 1)
              default_rr_eliminator <1>::process_diag_bhp (res, rr, dt);
          }
      }
    else
      {
        if (well->get_well_controller ()->is_rate ())
          {
            if (block_size == 3)
              default_rr_eliminator <3>::process_rate (res, rw, wr, dt);
            else if (block_size == 2)
              default_rr_eliminator <2>::process_rate (res, rw, wr, dt);
            else if (block_size == 1)
              default_rr_eliminator <1>::process_rate (res, rw, wr, dt);
          }
      }
  }

  void
  default_well::fill_rows (index_array_t &rows) const
  {
    detail::fill_rows <default_well, default_connection_t> (this, *rows);
  }
  void
  default_well::fill_jacobian (double dt, index_t block_size, const spv_long &rows, spv_long &cols, spv_float &values, stdv_long &markers) const
  {
    detail::fill_jacobian <default_well, default_connection_t> (this,
      dt, block_size, rows, cols, values, markers);
  }

  void
  default_well::fill_rhs (double dt, index_t n_phases, bool is_g, bool is_o, bool is_w, rhs_item_array_t &rhs_) const
  {
    item_t wefac = base_t::exploitation_factor_ > 0 ? base_t::exploitation_factor_ * dt : dt;

    typedef default_connection_iterator_impl <default_well, default_connection> iterator_t;
    iterator_t it (this, begin_iterator_tag), e (this, end_iterator_tag);
    rhs_item_t *rhs = &(*rhs_)[0];
    for (; !base_t::is_shut () && it != e; ++it)
      {
        const sp_connection_t &c = *it;
        if (!c->is_shut ())
          {
            int n_block = c->n_block ();
            int index   = n_block * n_phases;
            const shared_vector <rhs_item_t> &c_rate  = c->get_rate_value ();

            if (n_phases == 3)
              {
                rhs[index + p3_gas] += wefac * c_rate[p3_gas];
                rhs[index + p3_oil] += wefac * c_rate[p3_oil];
                rhs[index + p3_wat] += wefac * c_rate[p3_wat];
              }
            else if (n_phases == 2)
              {
                if (is_w)
                  {
                    rhs[index + p2ow_oil] += wefac * c_rate[p2ow_oil];
                    rhs[index + p2ow_wat] += wefac * c_rate[p2ow_wat];
                  }
                else if (is_g)
                  {
                    rhs[index + p2og_gas] += wefac * c_rate[p2og_gas];
                    rhs[index + p2og_oil] += wefac * c_rate[p2og_oil];
                  }
              }
            else
              {
                if (is_w)
                  rhs[index] += wefac * c_rate[0];
                else if (is_g)
                  rhs[index] += wefac * c_rate[0];
                else
                  {
                    BS_ASSERT (is_o);
                    rhs[index] += wefac * c_rate[0];
                  }
              }
          }
      }
  }

  template <bool is_prod>
  void
  default_well::calc_rate_and_derivs (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, BS_SP (jacobian) &jacobian)
  {
    if (calc_model->n_phases == 3)
      {
        calc_rate_and_derivs_concrete <true, true, true, is_prod> (calc_model, mesh, jacobian);
      }
    else if (calc_model->is_water () && calc_model->is_oil ())
      {
        calc_rate_and_derivs_concrete <true, false, true, is_prod> (calc_model, mesh, jacobian);
      }
    else if (calc_model->is_gas () && calc_model->is_oil ())
      {
        calc_rate_and_derivs_concrete <false, true, true, is_prod> (calc_model, mesh, jacobian);
      }
    else if (calc_model->is_water ())
      {
        calc_rate_and_derivs_concrete <true, false, false, is_prod> (calc_model, mesh, jacobian);
      }
    else if (calc_model->is_gas ())
      {
        calc_rate_and_derivs_concrete <false, true, false, is_prod> (calc_model, mesh, jacobian);
      }
    else
      {
        calc_rate_and_derivs_concrete <false, false, true, is_prod> (calc_model, mesh, jacobian);
      }
  }

  /**
   * \brief  Prints rr, pw, rw, wr values of well t
   * \param  t Instance of default_well
   * */
  template <typename T>
  void
  print (const T &t)
  {
    std::cout << "block[" << t->n_block () << "]:" << std::endl;
    //std::cout << "\tprod: " << std::endl;
    //std::cout << "\t\tw = " << t->rate_.prod.water << "\n\t\tg = " << t->rate_.prod.gas << "\n\t\to = " << t->rate_.prod.oil << std::endl;
    //std::cout << "\tinj: " << std::endl;
    //std::cout << "\t\tw = " << t->rate_.inj.water << "\n\t\tg = " << t->rate_.inj.gas << "\n\t\to = " << t->rate_.inj.oil << std::endl;
    //std::cout << "\tmob: " << std::endl;
    //std::cout << "\t\tw = " << t->mobility_value[0] << "\n\t\tg = " << t->mobility_value[1] << "\n\t\to = " << t->mobility_value[2] << std::endl;
    for (size_t i = 0, cnt = t->get_rr_value ().size (); i < cnt; ++i)
      {
        printf ("rr[%d] = %10.15f\n", i, t->get_rr_value ()[i]);
      }
    for (size_t i = 0, cnt = t->get_ps_value ().size (); i < cnt; ++i)
      {
        printf ("ps[%d] = %10.15f\n", i, t->get_ps_value ()[i]);
      }
    for (size_t i = 0, cnt = t->get_rw_value ().size (); i < cnt; ++i)
      {
        printf ("rw[%d] = %10.15f\n", i, t->get_rw_value ()[i]);
      }
    for (size_t i = 0, cnt = t->get_wr_value ().size (); i < cnt; ++i)
      {
        printf ("wr[%d] = %10.15f\n", i, t->get_wr_value ()[i]);
      }
  }
  /**
   * \brief  Prints production and injection rates of T
   * \param  t Instance of default_well or default_connection
   * */
  template <typename T>
  void
  print_rate (const T &t)
  {
    std::cout << "block[" << t->n_block () << "]:" << std::endl;
    std::cout << "\tprod: " << std::endl;
    std::cout << "\t\tw = " << t->rate_.prod.water << "\n\t\tg = " << t->rate_.prod.gas << "\n\t\to = " << t->rate_.prod.oil << std::endl;
    std::cout << "\tinj: " << std::endl;
    std::cout << "\t\tw = " << t->rate_.inj.water << "\n\t\tg = " << t->rate_.inj.gas << "\n\t\to = " << t->rate_.inj.oil << std::endl;
    //std::cout << "\tmob: " << std::endl;
    //std::cout << "\t\tw = " << t->mobility_value[0] << "\n\t\tg = " << t->mobility_value[1] << "\n\t\to = " << t->mobility_value[2] << std::endl;
  }

  template <bool is_w, bool is_g, bool is_o, bool is_prod>
  void
  default_well::calc_rate_and_derivs_concrete (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, BS_SP (jacobian) &jacobian)
  {

    calc_rate_and_derivs_t <is_w, is_g, is_o, is_prod> calc_ (calc_model, jacobian->get_sp_diagonal (), jacobian->get_sec_rhs (), this);

    //clear_data ();
    //for (size_t i = 0; i < 10000; ++i)
    //  {
    //    calc_.calc_rate (this);
    //    //base_t::well_controller_->calc_rate (calc_model, sp_this, jmatrix);
    //  }

    smart_ptr <well_t> sp_this (this);

    clear_data ();
    calc_.calc_rate (this);
    //base_t::well_controller_->calc_rate (calc_model, sp_this, jmatrix);
    //for (size_t i = 0, cnt = connection_list_.size (); i < cnt; ++i)
    //  {
    //    print_rate (connection_list_[i]);
    //  }

    if (base_t::well_controller_->check (sp_this))
      {
        base_t::calc_perf_bhp_->calculate (sp_this, calc_model, mesh);
        base_t::well_state_.is_work = check_connections_bhp (calc_model->pressure);
        if (!base_t::well_state_.is_work)
          {
            BOSOUT (section::wells, level::medium) << "[" << base_t::name_ << "] well does not work (2)" << bs_end;
            //init_approx_is_calc_ = false;
            return ;
          }

        clear_data ();
        calc_.calc_rate (this);
        //base_t::well_controller_->calc_rate (calc_model, sp_this, jmatrix);
        //for (size_t i = 0, cnt = connection_list_.size (); i < cnt; ++i)
        //  {
        //    print_rate (connection_list_[i]);
        //  }
      }

    calc_.calc_derivs (this, base_t::is_rate ());
    //base_t::well_controller_->calc_derivs (calc_model, sp_this, jmatrix);
    //for (size_t i = 0, cnt = connection_list_.size (); i < cnt; ++i)
    //  {
    //    print (connection_list_[i]);
    //  }
  }


  void
  default_well::process_internal (bool is_start, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, BS_SP (jacobian) &jacobian)
  {
    BS_ASSERT (!base_t::is_shut ()) (base_t::name ());

    smart_ptr <well_t> sp_this (this);

    if (is_start)
      {
        BS_ASSERT (base_t::calc_perf_density_);
        base_t::calc_perf_density_->calculate (sp_this, calc_model);
      }
    clear_data ();

    if (!base_t::init_approx_is_calc_)
      {
        BS_ASSERT (base_t::calc_rho_);
        BS_ASSERT (base_t::calc_well_pressure_);

        bool switched_to_bhp = base_t::calc_well_pressure_->calculate (sp_this, calc_model);
        if (switched_to_bhp)
          {
            base_t::calc_perf_density_->calculate (sp_this, calc_model);
          }

        base_t::init_approx_is_calc_ = true;
      }

    BS_ASSERT (base_t::calc_perf_bhp_);
    base_t::calc_perf_bhp_->calculate (sp_this, calc_model, mesh);
    base_t::well_state_.is_work = check_connections_bhp (calc_model->pressure);
    if (!base_t::well_state_.is_work)
      {
        BOSOUT (section::wells, level::medium) << "[" << base_t::name_ << "] well does not work" << bs_end;
        //init_approx_is_calc_ = false;
        return ;
      }

    if (base_t::well_controller_->is_production ())
      {
        calc_rate_and_derivs <true> (calc_model, mesh, jacobian);
      }
    else
      {
        calc_rate_and_derivs <false> (calc_model, mesh, jacobian);
      }
  }

  void
  default_well::process_impl (bool is_start, double /*dt*/, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, BS_SP (jacobian) &jacobian)
  {
    process_internal (is_start, calc_model, mesh, jacobian);

    if (base_t::is_rate ())
      {
        BOSOUT (section::wells, level::low) << "[" << this->name_ << "] start newton iterations" << bs_end;

        static index_t bw_iter_max = 20;
        static item_t bw_eps = 1.0e-7f;
        //item_t prev_bw_value = bw_value;

        for (index_t bw_iter = 0; bw_iter < bw_iter_max; ++bw_iter)
          {
            if (fabs (bw_value) <= bw_eps)
              {
#ifdef _DEBUG
                BOSOUT (section::wells, level::debug) << "[" << base_t::name_ << "] bw_value (" << bw_value << ") smaller than eps" << bs_end;
#endif
                break;
              }

            ww_value += bw_eps;
            //prev_bw_value = bw_value;
            item_t bw_ww_dx = bw_value / ww_value;
#ifdef _DEBUG
            BOSOUT (section::wells, level::debug) << "[" << this->name_ << "] bw: " << bw_value << bs_end;
            BOSOUT (section::wells, level::debug) << "[" << this->name_ << "] ww: " << ww_value << bs_end;
#endif
            base_t::bhp_ += bw_ww_dx;

            process_internal (false, calc_model, mesh, jacobian);
            if (!base_t::is_rate ())
              break;
          }

        BOSOUT (section::wells, level::low) << "[" << this->name_ << "] end newton iterations" << bs_end;
      }
  }

  void
  default_well::clear_data ()
  {
    typedef default_connection_iterator_impl <default_well, default_connection> iterator_t;
    iterator_t it (iterator_t (this, begin_iterator_tag)), 
               e (iterator_t (this, end_iterator_tag));
    for (; it != e; ++it)
      {
        it->clear_data ();
      }

    ww_value = 0;
    bw_value = 0;

    base_t::clear_data ();
  }

  bool
  default_well::check_shut ()
  {
    open_connections_count_ = 0;
    if (base_t::is_shut ())
      {
        return true;
      }

    typedef default_connection_iterator_impl <default_well, default_connection> iterator_t;
    iterator_t it (this, begin_iterator_tag), e (this, end_iterator_tag);
    for (; it != e; ++it)
      {
        if (!it->is_shut ())
          {
            open_connections_count_++;
          }
      }

    return false;
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (default_well);
  BLUE_SKY_TYPE_STD_COPY (default_well);
  BLUE_SKY_TYPE_IMPL (default_well, well, "default_well", "default_well", "default_well");
  //////////////////////////////////////////////////////////////////////////

  bool
  default_well_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, default_connection::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, default_well::bs_type ()); BS_ASSERT (res);

    return res;
  }

} // namespace wells
} // namespace blue_sky


