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
    connection::connection (blue_sky::bs_type_ctor_param param /* = NULL */)
    : n_block_ (-1)
    {
    }
    /**
     * \brief  copy-ctor for connection
     * \param  c Instance of connection to be copied
     * */
    connection::connection (const connection_t &c)
    : bs_refcounter (c), objbase (c)
    {
      *this = c;
    }

    shared_vector <t_double>
    connection::get_rw_value ()
    {
      return shared_array <item_t> (0, 0);
    }
    shared_vector <t_double>
    connection::get_wr_value ()
    {
      return shared_array <item_t> (0, 0);
    }
    shared_vector <t_double>
    connection::get_rr_value ()
    {
      return shared_array <item_t> (0, 0);
    }
    shared_vector <t_double>
    connection::get_ps_value ()
    {
      return shared_array <item_t> (0, 0);
    }
    shared_vector <t_double>
    connection::get_rate_value ()
    {
      return shared_array <rhs_item_t> (0, 0);
    }

    connection::item_t
    connection::get_cur_bhp () const
      {
        return cur_bhp;
      }
    connection::item_t
    connection::get_density () const
      {
        return density;
      }
    connection::item_t
    connection::get_connection_depth () const
      {
        return connection_depth;
      }
    connection::item_t
    connection::get_bulkp () const
      {
        return bulkp;
      }

    //////////////////////////////////////////////////////////////////////////

    void
    connection::compute_factors (const physical_constants &internal_constants,
        const sp_params_t &params,
        const sp_mesh_iface_t &mesh,
        const stdv_float &perm,
        const stdv_float &ntg,
        bool ro_calc_flag)
    {
      BS_ASSERT (n_block_ >= 0) (n_block_);
      if (n_block_ < 0)
        throw bs_exception ("connection::compute_factors_by_peaceman_model", "Invalid block number");

      using namespace wells::compute_factors;

      if (HORIZ_WELL_MODEL_PEACEMAN || dir_ == direction_z /* || !wfrictn_keyword_active*/)
        {
          peaceman_model::compute (*this, internal_constants, params, mesh, perm, ntg, ro_calc_flag);
        }
      //else
      //  baby_odeh_model::compute (*this, internal_contstants, d1, d2, d3, perm1, perm2, perm3, ntg, ro_calc_flag);
    }

    void
    connection::mul_perm_mult (item_t mult)
    {
      //BOSOUT (section::wells, level::debug)
      //  << "WPIMULT: " << (*it)->n_block ()
      //  << "  before: " << (*it)->mult ()
      //  << "  after:  " << (*it)->mult () * perm_mult
      //  << bs_end;

      mult_ *= mult;
    }
    void
    connection::set_half_length (item_t half_length)
    {
      fracture_half_length_ = half_length;
    }
    void
    connection::set_theta (item_t angle)
    {
      fracture_angle_ = angle;
    }
    void
    connection::set_skin (item_t skin)
    {
      skin_ = skin;
    }
    void
    connection::set_status (connection_status_type connection_status)
    {
      status_ = connection_status;
    }
    void
    connection::set_factor (item_t factor)
    {
      fact_ = factor;
    }
    void
    connection::set_diameter (item_t diameter)
    {
      diam_ = diameter;
    }
    void
    connection::set_Kh (item_t kh)
    {
      kh_ = kh;
    }
    void
    connection::set_direction (connection_direction_type direction)
    {
      dir_ =  direction;
    }
    void
    connection::set_coord (index_t i, index_t j, index_t k, index_t nb)
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
    void
    connection::set_seg_number (const index_t &seg)
    {
      iseg_ = seg;
    }

    t_long
    connection::n_block () const
    {
      return n_block_;
    }
    t_long
    connection::i_coord () const
    {
      return i_coord_;
    }
    t_long
    connection::j_coord () const
    {
      return j_coord_;
    }
    t_long
    connection::k_coord () const
    {
      return k_coord_;
    }
    t_double
    connection::mult () const
    {
      return mult_;
    }

    t_double
    connection::get_head_term () const
    {
      return head_term;
    }
    t_long
    connection::get_seg_number () const
    {
      return iseg_;
    }

    t_double
    connection::get_fact () const
      {
        return fact_;
      }

    void
    connection::set_cur_bhp (item_t bhp)
    {
      cur_bhp = bhp;
    }
    void
    connection::set_bulkp (item_t b)
    {
      bulkp = b;
    }
    void
    connection::set_rate (item_t rate)
    {
      BS_ASSERT (false && "NOT IMPL YET");
    }
    void
    connection::set_head_term (item_t term)
    {
      head_term = term;
    }
    void
    connection::set_mult (item_t mult)
    {
      mult_ = mult;
    }

    bool
    connection::is_shut () const
      {
        return status_ == connection_shut;
      }

    void
    connection::clear_data ()
    {
      rate_     = 0;
      rate_rc_  = 0;
    }

    void
    connection::set_connection_depth (const sp_mesh_iface_t &mesh)
    {
      BS_ASSERT (i_coord_ != -1 && j_coord_ != -1 && k_coord_ != -1) (i_coord_) (j_coord_) (k_coord_);
      BS_ASSERT (n_block_ != -1) (n_block_);
      connection_depth = mesh->get_element_depth (n_block_);
      //BOSOUT (section::wells, level::debug) << "[xxx : " << n_block_ << "] set_connection_depth: " << (index_t)i_coord_ << " " << (index_t)j_coord_ << " " << (index_t)k_coord_
      //  << " = " << connection_depth
      //  << bs_end;
    }


    BLUE_SKY_TYPE_STD_CREATE (connection)
    BLUE_SKY_TYPE_STD_COPY (connection)
    BLUE_SKY_TYPE_IMPL (connection, objbase, "well_connection_seq", "BOS_Core well connection class", "BOS_Core well connection class")

  } // namespace wells
} // namespace blue_sky
