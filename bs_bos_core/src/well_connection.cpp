/**
 *       \file  well_connection.cpp
 *      \brief  Implementation of base class for well
 *              perforations (connections)
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  15.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "calc_well.h"
#include "wells_compute_connection_factors.h"
#include "well_connection.h"
#include "calc_model.h"

#define HORIZ_WELL_MODEL_PEACEMAN 1

namespace blue_sky
  {
  namespace wells
    {

    ///////////////////////////////////////////////////////////////////////////
    connection_direction_type
    connection_direction_cast (const std::string &str)
    {
      if (str == "X" || str == "x")
        return direction_x;
      else if (str == "Y" || str == "y")
        return direction_y;
      else if (str == "Z" || str == "z")
        return direction_z;
      else if (str == "")
        return direction_z;
      else
        {
          BS_ASSERT (false && "Unsupported value of direction") (str);
          throw bs_exception ("connection_direction_cast", "Unsupported value of direction");
        }
    }
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    connection_status_type
    connection_status_cast (const std::string &str)
    {
      if (str == "OPEN")
        return connection_open;
      else if (str == "SHUT")
        return connection_shut;
      else if (str == "")
        return connection_shut;
      else
        {
          BS_ASSERT (false && "Unsupported value of connection status") (str);
          throw bs_exception ("connection_status_cast", "Unsupported value of connection");
        }
    }
    ///////////////////////////////////////////////////////////////////////////

    /**
     * \brief  'default' ctor for connection
     * \param  param Additional params for ctor
     * */
    template <typename strategy_t>
    connection<strategy_t>::connection (blue_sky::bs_type_ctor_param param /* = NULL */)
    : n_block_ (-1)
    {
    }
    /**
     * \brief  copy-ctor for connection
     * \param  c Instance of connection to be copied
     * */
    template <typename strategy_t>
    connection<strategy_t>::connection (const connection_t &c)
    : bs_refcounter (c), objbase (c)
    {
      *this = c;
    }

    template <typename strategy_t>
    array_ext <typename strategy_t::item_t>
    connection<strategy_t>::get_rw_value ()
    {
      return array_ext <item_t> (0, 0);
    }
    template <typename strategy_t>
    array_ext <typename strategy_t::item_t>
    connection<strategy_t>::get_wr_value ()
    {
      return array_ext <item_t> (0, 0);
    }
    template <typename strategy_t>
    array_ext <typename strategy_t::item_t>
    connection<strategy_t>::get_rr_value ()
    {
      return array_ext <item_t> (0, 0);
    }
    template <typename strategy_t>
    array_ext <typename strategy_t::item_t>
    connection<strategy_t>::get_ps_value ()
    {
      return array_ext <item_t> (0, 0);
    }
    template <typename strategy_t>
    array_ext <typename strategy_t::rhs_item_t>
    connection<strategy_t>::get_rate_value ()
    {
      return array_ext <rhs_item_t> (0, 0);
    }

    template <typename strategy_t>
    typename connection<strategy_t>::item_t
    connection<strategy_t>::get_cur_bhp () const
      {
        return cur_bhp;
      }
    template <typename strategy_t>
    typename connection<strategy_t>::item_t
    connection<strategy_t>::get_density () const
      {
        return density;
      }
    template <typename strategy_t>
    typename connection<strategy_t>::item_t
    connection<strategy_t>::get_connection_depth () const
      {
        return connection_depth;
      }
    template <typename strategy_t>
    typename connection<strategy_t>::item_t
    connection<strategy_t>::get_bulkp () const
      {
        return bulkp;
      }

    //////////////////////////////////////////////////////////////////////////

    template <typename strategy_t>
    void
    connection<strategy_t>::compute_factors (const physical_constants &internal_contstants,
        const sp_params_t &params,
        const sp_mesh_iface_t &mesh,
        const item_array_t &perm,
        const item_array_t &ntg,
        bool ro_calc_flag)
    {
      BS_ASSERT (n_block_ >= 0) (n_block_);
      if (n_block_ < 0)
        throw bs_exception ("connection::compute_factors_by_peaceman_model", "Invalid block number");

      using namespace wells::compute_factors;

      if (HORIZ_WELL_MODEL_PEACEMAN || dir_ == direction_z /* || !wfrictn_keyword_active*/)
        {
          peaceman_model<strategy_t>::compute (*this, internal_contstants, params, mesh, perm, ntg, ro_calc_flag);
        }
      //else
      //  baby_odeh_model::compute (*this, internal_contstants, d1, d2, d3, perm1, perm2, perm3, ntg, ro_calc_flag);
    }

    template <typename strategy_t>
    void
    connection<strategy_t>::mul_perm_mult (item_t mult)
    {
      //BOSOUT (section::wells, level::debug)
      //  << "WPIMULT: " << (*it)->n_block () 
      //  << "  before: " << (*it)->mult ()
      //  << "  after:  " << (*it)->mult () * perm_mult
      //  << bs_end; 

      mult_ *= mult;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_half_length (item_t half_length)
    {
      fracture_half_length_ = half_length;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_theta (item_t angle)
    {
      fracture_angle_ = angle;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_skin (item_t skin)
    {
      skin_ = skin;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_status (connection_status_type connection_status)
    {
      status_ = connection_status;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_factor (item_t factor)
    {
      fact_ = factor;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_diameter (item_t diameter)
    {
      diam_ = diameter;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_Kh (item_t kh)
    {
      kh_ = kh;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_direction (connection_direction_type direction)
    {
      dir_ =  direction;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_coord (index_t i, index_t j, index_t k, index_t nb)
    {
      BS_ASSERT (i >= 0) (i);
      BS_ASSERT (j >= 0) (j);
      BS_ASSERT (k >= 0) (k);
      BS_ASSERT (nb >= 0) (nb);

      i_coord_ = i;
      j_coord_ = j;
      k_coord_ = k;
      n_block_ = nb;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_seg_number (const index_t &seg)
    {
      iseg_ = seg;
    }

    template <typename strategy_t>
    typename connection<strategy_t>::index_t
    connection<strategy_t>::n_block () const
    {
      return n_block_;
    }
    template <typename strategy_t>
    typename connection<strategy_t>::index_t
    connection<strategy_t>::i_coord () const
    {
      return i_coord_;
    }
    template <typename strategy_t>
    typename connection<strategy_t>::index_t
    connection<strategy_t>::j_coord () const
    {
      return j_coord_;
    }
    template <typename strategy_t>
    typename connection<strategy_t>::index_t
    connection<strategy_t>::k_coord () const
    {
      return k_coord_;
    }
    template <typename strategy_t>
    typename connection<strategy_t>::item_t
    connection<strategy_t>::mult () const
    {
      return mult_;
    }

    template <typename strategy_t>
    typename connection <strategy_t>::item_t
    connection <strategy_t>::get_head_term () const
    {
      return head_term;
    }
    template <typename strategy_t>
    typename connection <strategy_t>::index_t
    connection <strategy_t>::get_seg_number () const
    {
      return iseg_;
    }

    template <typename strategy_t>
    typename strategy_t::item_t
    connection<strategy_t>::get_fact () const
      {
        return fact_;
      }

    template <typename strategy_t>
    void
    connection<strategy_t>::set_cur_bhp (item_t bhp)
    {
      cur_bhp = bhp;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_bulkp (item_t b)
    {
      bulkp = b;
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_rate (item_t rate)
    {
      BS_ASSERT (false && "NOT IMPL YET");
    }
    template <typename strategy_t>
    void
    connection<strategy_t>::set_head_term (item_t term)
    {
      head_term = term;
    }
    template <typename strategy_t>
    void
    connection <strategy_t>::set_mult (item_t mult)
    {
      mult_ = mult;
    }

    template <typename strategy_t>
    bool
    connection <strategy_t>::is_shut () const
      {
        return status_ == connection_shut;
      }

    template <typename strategy_t>
    void
    connection<strategy_t>::clear_data ()
    {
      rate_     = 0;
      rate_rc_  = 0;
    }

    template <typename strategy_t>
    void
    connection<strategy_t>::set_connection_depth (const sp_mesh_iface_t &mesh)
    {
      BS_ASSERT (i_coord_ != -1 && j_coord_ != -1 && k_coord_ != -1) (i_coord_) (j_coord_) (k_coord_);
      BS_ASSERT (n_block_ != -1) (n_block_);
      connection_depth = mesh->get_element_depth (n_block_);
      //BOSOUT (section::wells, level::debug) << "[xxx : " << n_block_ << "] set_connection_depth: " << (index_t)i_coord_ << " " << (index_t)j_coord_ << " " << (index_t)k_coord_
      //  << " = " << connection_depth
      //  << bs_end;
    }


    BLUE_SKY_TYPE_STD_CREATE_T_DEF (connection, (class))
    BLUE_SKY_TYPE_STD_COPY_T_DEF (connection, (class))
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (connection <base_strategy_fi>), 1, (objbase), "well_connection_seq_fi", "BOS_Core well connection class", "BOS_Core well connection class", false)
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (connection <base_strategy_di>), 1, (objbase), "well_connection_seq_di", "BOS_Core well connection class", "BOS_Core well connection class", false)
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (connection <base_strategy_mixi>), 1, (objbase), "well_connection_seq_mixi", "BOS_Core well connection class", "BOS_Core well connection class", false)

  } // namespace wells
} // namespace blue_sky
