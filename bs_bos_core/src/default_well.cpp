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
#include "well_rate_control_deriv.h"
#include "default_well_calc_rate.h"

namespace blue_sky {
namespace wells {

  /**
   * \brief  'default' ctor for default_well
   * \param  param Additional params for ctor
   * */
  template <typename strategy_t>
  default_well <strategy_t>::default_well (bs_type_ctor_param param /* = NULL */)
  : base_t (param)
  , ww_value (0)
  , bw_value (0)
  {
  }

  /**
   * \brief  ctor for default_well
   * \param  well_name Name of well
   * \todo   Obsolete 
   * */
  template <typename strategy_t>
  default_well <strategy_t>::default_well (const std::string &well_name)
  : base_t (well_name)
  , ww_value (0)
  , bw_value (0)
  {

  }

  /**
   * \brief  copy-ctor for default_well
   * \param  w Instance of default_well to be copied
   * */
  template <typename strategy_t>
  default_well <strategy_t>::default_well (const default_well &w)
        : bs_refcounter (), well < strategy_t> ()
  {
    *this = w;
  }

  namespace detail {

    /**
     * \class sort_connection_list
     * \brief Sorts connection list in ASC order
     * */
    template <typename connection_t>
    struct sort_connection_list : std::binary_function <const connection_t &, const connection_t &, bool>
    {
      /**
       * \brief  Sorts connection list in ASC order
       * \param  lhs
       * \param  rhs
       * \return True if n_block of lhs < of n_block of rhs
       * */
      bool
      operator () (const connection_t &lhs, const connection_t &rhs) const
      {
        return lhs->n_block () < rhs->n_block ();
      }
    };
  } // namespace wells

  template <typename strategy_t>
  typename default_well <strategy_t>::sp_connection_t
  default_well <strategy_t>::add_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block)
  {
    if (n_block < 0)
      bs_throw_exception ("Invalid connection n_block value");

    sp_default_connection_t connection = BS_KERNEL.create_object (default_connection_t::bs_type (), true);
    if (!connection)
        bs_throw_exception ("Can't create connection");

    connection->set_coord (i_coord, j_coord, k_coord, n_block);
    connection_list_.push_back (connection);
    std::sort (connection_list_.begin (), connection_list_.end (), detail::sort_connection_list <sp_default_connection_t> ());
    connection_map_.insert (std::make_pair (connection->n_block (), (index_t)connection_list_.size () - 1));

    return connection;
  }

  template <typename strategy_t>
  typename default_well <strategy_t>::sp_connection_t
  default_well <strategy_t>::get_connection (index_t idx) const
  {
    BS_ASSERT (idx < connection_list_.size ()) (base_t::name ()) (idx) (connection_list_.size ());
    return connection_list_[idx];
  }
  template <typename strategy_t>
  typename default_well <strategy_t>::sp_connection_t
  default_well <strategy_t>::get_connection_map (index_t n_block) const
  {
    for (size_t i = 0, cnt = connection_list_.size (); i < cnt; ++i)
      {
        if (connection_list_[i]->n_block () == n_block)
          return connection_list_[i];
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

  template <typename strategy_t>
  size_t
  default_well <strategy_t>::get_connections_count () const
  {
    return connection_list_.size ();
  }


  template <typename strategy_t>
  void
  default_well <strategy_t>::restore_solution (double /*dt*/, const item_array_t &p_sol, 
                                               const item_array_t & /*s_sol*/, index_t block_size)
  {
    BS_ASSERT (get_connections_count ()) (base_t::name ());

    typedef array_ext <item_t>        item_wr_block_t;
    typedef array_ext <const item_t>  item_xr_block_t;

    // TODO: may be we can not check is_work flag
    if (!base_t::well_state_.is_work)
      return ;

    if (!base_t::well_controller_->is_rate ())
      return ;

    item_t xw = 0;
    item_t ww_bw = (1.0 / ww_value) * bw_value;

    for (index_t j = 0, jcnt = (index_t)base_t::open_connections_.size (); j < jcnt; ++j)
      {
        index_t con_index = base_t::open_connections_[j];
        BS_ASSERT (connection_list_[con_index]->is_shut () == false) (j);

        const sp_connection_t &c  = connection_list_[con_index];
        const item_wr_block_t &wr = c->get_wr_value ();
        const item_xr_block_t &xr = item_xr_block_t (&p_sol[c->n_block () * block_size], block_size);

        item_t wr_xr = 0;
        for (index_t j = 0; j < block_size; ++j)
          {
            wr_xr += wr[j] * xr[j];
          }

        xw += ww_bw - wr_xr;
      }

    BOSOUT (section::wells, level::low) << "[" << base_t::name_ << boost::format ("] Restore solution: change bhp from %f to %f (bw: %f, ww: %f)") % base_t::bhp_ % (base_t::bhp_ + xw) % bw_value % ww_value << bs_end;
    base_t::bhp_ += xw;
  }
  template <typename strategy_t>
  void
  default_well <strategy_t>::eliminate (rhs_item_t *res, index_t rw_index, index_t wr_index, double dt, index_t block_size) const
  {
    BS_ASSERT (!base_t::is_shut ()) (base_t::name ());
    BS_ASSERT (get_connections_count ()) (base_t::name ());
    BS_ASSERT (res);

    if (block_size > 3 && block_size < 1)
      bs_throw_exception ("block_size invalid");

    const sp_connection_t &rw_con = connection_list_[rw_index];
    const sp_connection_t &wr_con = connection_list_[wr_index];

    const array_ext <item_t> &rw = rw_con->get_rw_value ();
    const array_ext <item_t> &wr = wr_con->get_wr_value ();
    const array_ext <item_t> &rr = wr_con->get_rr_value ();

    BS_ASSERT (rw.size () == (size_t)block_size) (rw.size ()) (block_size);
    BS_ASSERT (wr.size () == (size_t)block_size) (wr.size ()) (block_size);
    BS_ASSERT (rr.size () == (size_t)block_size * block_size) (rr.size ()) (block_size);

    if (rw_index == wr_index)
      {
        if (well_t::well_controller_->is_rate ())
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
        if (well_t::well_controller_->is_rate ())
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

  template <typename strategy_t>
  template <bool is_prod>
  void
  default_well <strategy_t>::calc_rate_and_derivs (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix)
  {
    if (calc_model->n_phases == 3)
      {
        calc_rate_and_derivs_concrete <true, true, true, is_prod> (calc_model, mesh, jmatrix);
      }
    else if (calc_model->is_water () && calc_model->is_oil ())
      {
        calc_rate_and_derivs_concrete <true, false, true, is_prod> (calc_model, mesh, jmatrix);
      }
    else if (calc_model->is_gas () && calc_model->is_oil ())
      {
        calc_rate_and_derivs_concrete <false, true, true, is_prod> (calc_model, mesh, jmatrix);
      }
    else if (calc_model->is_water ())
      {
        calc_rate_and_derivs_concrete <true, false, false, is_prod> (calc_model, mesh, jmatrix);
      }
    else if (calc_model->is_gas ())
      {
        calc_rate_and_derivs_concrete <false, true, false, is_prod> (calc_model, mesh, jmatrix);
      }
    else
      {
        calc_rate_and_derivs_concrete <false, false, true, is_prod> (calc_model, mesh, jmatrix);
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

  template <typename strategy_t>
  template <bool is_w, bool is_g, bool is_o, bool is_prod>
  void
  default_well <strategy_t>::calc_rate_and_derivs_concrete (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix)
  {

    calc_rate_and_derivs_t <strategy_t, is_w, is_g, is_o, is_prod> calc_ (calc_model, jmatrix, this, this->connection_list_);

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


  template <typename strategy_t>
  void
  default_well <strategy_t>::process_internal (bool is_start, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix)
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
        calc_rate_and_derivs <true> (calc_model, mesh, jmatrix);
      }
    else
      {
        calc_rate_and_derivs <false> (calc_model, mesh, jmatrix);
      }
  }

  template <typename strategy_t>
  void
  default_well <strategy_t>::process (bool is_start, double /*dt*/, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix)
  {
    if (base_t::is_shut ())
      {
        BOSOUT (section::wells, level::debug) << "[" << base_t::name () << "] is_shut" << bs_end;
        return ;
      }

#ifdef _DEBUG
    if (!base_t::is_shut () &&  connection_list_.empty ())
      {
        bs_throw_exception (boost::format ("[%s]: not shut but connection list is empty") % base_t::name ());
      }
#endif

    process_internal (is_start, calc_model, mesh, jmatrix);

    if (base_t::is_rate ())
      {
        BOSOUT (section::wells, level::low) << "[" << this->name_ << "] start newton iterations" << bs_end;

        static index_t bw_iter_max = 20;
        static item_t bw_eps = 1.0e-7f;
        item_t prev_bw_value = bw_value;

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
            prev_bw_value = bw_value;
            item_t bw_ww_dx = bw_value / ww_value;
#ifdef _DEBUG
            BOSOUT (section::wells, level::debug) << "[" << this->name_ << "] bw: " << bw_value << bs_end;
            BOSOUT (section::wells, level::debug) << "[" << this->name_ << "] ww: " << ww_value << bs_end;
#endif
            base_t::bhp_ += bw_ww_dx;

            process_internal (false, calc_model, mesh, jmatrix);
            if (!base_t::is_rate ())
              break;
          }

        BOSOUT (section::wells, level::low) << "[" << this->name_ << "] end newton iterations" << bs_end;
      }
  }

  template <typename strategy_t>
  void
  default_well <strategy_t>::clear_data ()
  {
    for (size_t i = 0, cnt = connection_list_.size (); i < cnt; ++i)
      {
        connection_list_ [i]->clear_data ();
      }

    ww_value = 0;
    bw_value = 0;

    base_t::clear_data ();
  }

  template <typename strategy_t>
  array_ext <typename default_well <strategy_t>::item_t>
  default_well <strategy_t>::get_ww_value ()
  {
    return array_ext <item_t> (&ww_value, 1);
  }
  template <typename strategy_t>
  array_ext <typename default_well <strategy_t>::item_t>
  default_well <strategy_t>::get_bw_value ()
  {
    return array_ext <item_t> (&bw_value, 1);
  }
  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (default_well, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (default_well, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_well<base_strategy_fi>), 1, (well<base_strategy_fi>), "default_well_seq_fi", "default_well_seq_fi", "default_well_seq_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_well<base_strategy_di>), 1, (well<base_strategy_di>), "default_well_seq_di", "default_well_seq_di", "default_well_seq_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_well<base_strategy_mixi>), 1, (well<base_strategy_mixi>), "default_well_seq_mixi", "default_well_seq_mixi", "default_well_seq_mixi", false);
  //////////////////////////////////////////////////////////////////////////

  bool
  default_well_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, default_connection<base_strategy_fi>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, default_connection<base_strategy_di>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, default_connection<base_strategy_mixi>::bs_type ()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, default_well<base_strategy_fi>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, default_well<base_strategy_di>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, default_well<base_strategy_mixi>::bs_type ()); BS_ASSERT (res);

    return res;
  }

} // namespace wells
} // namespace blue_sky


