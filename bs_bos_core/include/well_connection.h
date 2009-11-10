/**
 * \file well_connection.h
 * \brief well connection class declaration
 * \author Sergey Miryanov
 * \date 06.08.2008
 * */
#ifndef BS_WELL_CONNECTION_H_
#define BS_WELL_CONNECTION_H_

#include "fi_params.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "convert_units.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "well_type_helper.h"
#include "wells_compute_connection_factors.h"
#include "array_ext.h"

#include "well_controller.h"


namespace blue_sky
  {

  template <typename strategy_t>
  class well;

  template <typename strategy_t>
  class calc_model;

  namespace wells
    {

    ///////////////////////////////////////////////////////////////////////////
    enum connection_type
    {
      CONNECTION_USUAL,
      CONNECTION_GRP_PRIMARY,
      CONNECTION_GRP_SECONDARY,

      CONNECTION_TOTAL,
    };
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    enum connection_direction_type
    {
      direction_x,
      direction_y,
      direction_z,

      direction_total,
    };

    connection_direction_type
    connection_direction_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    enum connection_status_type
    {
      connection_open,
      connection_shut,

      connection_total,
    };

    connection_status_type
    connection_status_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    template <typename strategy_t>
    class BS_API_PLUGIN connection : public objbase
      {
      public:

        typedef connection <strategy_t>                   this_t;
        typedef smart_ptr <this_t, true>                  sp_this_t;

        typedef typename strategy_t::item_array_t         item_array_t;
        typedef typename strategy_t::item_t               item_t;
        typedef typename strategy_t::rhs_item_t           rhs_item_t;
        typedef typename strategy_t::index_t              index_t;

        typedef rate_data <strategy_t>                    rate_data_t;
        typedef typename rate_data_t::rate_data_inner     rate_data_inner_t;

        typedef fi_params                                 params_t;
        typedef smart_ptr <params_t>                      sp_params_t;

        typedef rs_mesh_iface <strategy_t>                mesh_iface_t;
        typedef smart_ptr <mesh_iface_t, true>            sp_mesh_iface_t;

        typedef connection <strategy_t>                   connection_t;
        typedef calc_model <strategy_t>                   calc_model_t;
        typedef well <strategy_t>                         well_t;

      public:

        void compute_factors (const physical_constants &internal_contstants,
                              const sp_params_t &params,
                              const sp_mesh_iface_t &mesh,
                              const item_array_t &perm,
                              const item_array_t &ntg,
                              bool ro_calc_flag = false);

        void mul_perm_mult (item_t mult);

        void set_half_length (item_t half_length);
        void set_theta (item_t theta);
        void set_skin (item_t skin);
        void set_status (connection_status_type connection_status);
        void set_factor (item_t factor);
        void set_diameter (item_t diameter);
        void set_Kh (item_t kh);
        void set_direction (connection_direction_type direction);
        void set_coord (index_t i, index_t j, index_t k, index_t n_block);
        void set_connection_depth (const sp_mesh_iface_t &mesh);

        index_t n_block () const;
        item_t get_fact () const;

        void set_bulkp (item_t bulkp);
        void set_rate (item_t rate);
        void set_head_term (item_t head_term);
        void set_cur_bhp (item_t cur_bhp);
        void set_mult (item_t mult);
        void set_seg_number (const index_t &seg);

        bool is_shut () const;

        connection_status_type get_status () const
          {
            return status_;
          }

        connection_direction_type get_dir () const
          {
            return dir_;
          }

        item_t                get_cur_bhp             () const;
        item_t                get_density             () const;
        item_t                get_connection_depth    () const;
        item_t                get_bulkp               () const;

        index_t               i_coord () const;
        index_t               j_coord () const;
        index_t               k_coord () const;

        item_t                mult () const;

        item_t                get_head_term () const;
        index_t               get_seg_number () const;

        virtual void clear_data ();
        virtual array_ext <item_t> get_rw_value   ();
        virtual array_ext <item_t> get_wr_value   ();
        virtual array_ext <item_t> get_rr_value   ();
        virtual array_ext <item_t> get_ps_value   ();
        virtual array_ext <rhs_item_t> get_rate_value ();

      public:
        auto_value <item_t>           head_term;
        auto_value <item_t>           cur_bhp;
        auto_value <item_t>           connection_depth;
        auto_value <item_t>           density;
        auto_value <item_t>						bulkp;

        rate_data_t                   rate_;
        rate_data_t                   rate_rc_;

        const rate_data_t &
        rate () const
        {
          return rate_;
        }

        const rate_data_inner_t &
        rate_prod () const
        {
          return rate_.prod;
        }

        const rate_data_inner_t &
        rate_inj () const
        {
          return rate_.inj;
        }

      private:

        template <typename strategy_t_x>
        friend struct wells::compute_factors::peaceman_model;

        template <typename strategy_t_x>
        friend struct wells::compute_factors::baby_odeh_model;

        auto_value <index_t, -1>			i_coord_;											//!< i coord of connection
        auto_value <index_t, -1>			j_coord_;											//!< j coord of connection
        auto_value <index_t, -1>			k_coord_;											//!< k coord of connection

        auto_value <connection_status_type, connection_open>
        status_;											//!< connection status

        // for compute factors. in the future this vars should be evolved to separate class
        auto_value <item_t>           diam_;                        //!< Well diameter
        auto_value <item_t>           kh_;                          //!< Effective Kh
        auto_value <item_t>           skin_;                        //!< Well skin factor
        auto_value <item_t, 1>        mult_;                        //!<
        auto_value <item_t, -1>       fact_;                        //!< connection factor
        auto_value <item_t, -1>       R0_;                          //!< R0

        auto_value <item_t>           fracture_half_length_;        //!< half length of fracture
        auto_value <item_t>           fracture_angle_;              //!< angle of fracture in grid block to OX direction
        auto_value <item_t>           fracture_perm_;               //!< coefficient of "permeability" on fracture

        auto_value <index_t, -1>      n_block_;                     //!< number of block in grid (n_block <-> (i, j, k))

        auto_value <connection_direction_type, direction_x>
        dir_;                         //!< direction

        auto_value <connection_type, CONNECTION_USUAL>
        connection_type_;             //!< connection type

      protected:
        auto_value <index_t, 0>       iseg_;
      public:
        BLUE_SKY_TYPE_DECL_T (connection <strategy_t>);
      };


  } // namespace wells
} // namespace blue_sky


#endif  // #ifndef BS_WELL_CONNECTION_H_
