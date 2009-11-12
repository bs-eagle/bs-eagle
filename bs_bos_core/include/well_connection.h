/**
 *       \file  well_connection.h
 *      \brief  Base class for well perforations (well connections)
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  06.08.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
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
    /**
     * \enum  connection_type
     * \brief Type of perforation (connection)
     * */
    enum connection_type
    {
      CONNECTION_USUAL,             //!< Usual perforation
      CONNECTION_GRP_PRIMARY,       //!< Primary GRP perforation
      CONNECTION_GRP_SECONDARY,     //!< Secondary GRP perforation

      CONNECTION_TOTAL,
    };
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    /**
     * \enum  connection_deriction_type
     * \brief Direction in which the well penetrates the grid block
     * */
    enum connection_direction_type
    {
      direction_x,                  //!< Along X axe
      direction_y,                  //!< Along Y axe
      direction_z,                  //!< Along Z axe

      direction_total,
    };

    /**
     * \brief  Converts string value to connection_direction_type
     * \param  str String value to convert
     * \return Throws exception if value is an invalid,
     *         direction_z if value is an empty otherwise
     *         element of connection_direction_type
     * */
    connection_direction_type
    connection_direction_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    /**
     * \enum  connection_status_type
     * \brief Perforation (connection) status
     * */
    enum connection_status_type
    {
      connection_open,            //!< Is perforation (connection) is open
      connection_shut,            //!< Is perforation (connection) is shut

      connection_total,
    };

    /**
     * \brief  Converts string value to connection_status_type
     * \param  str String value to convert
     * \return Throws exception if value is an invalid,
     *         connection_shut if value is an empty otherwise
     *         element of connection_status_type
     * */
    connection_status_type
    connection_status_cast (const std::string &str);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    /**
     * \class connection
     * \brief Base class for well perforations (well connections)
     * */
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

        /**
         * \brief  Computes perforation (connection) factor
         * \param  internal_constants
         * \param  params
         * \param  mesh
         * \param  perm
         * \param  ntg
         * \param  ro_calc_flag
         * */
        void 
        compute_factors (const physical_constants &internal_contstants,
                         const sp_params_t &params,
                         const sp_mesh_iface_t &mesh,
                         const item_array_t &perm,
                         const item_array_t &ntg,
                         bool ro_calc_flag = false);

        /**
         * \brief  Multiplies mult_ with mult
         * \param  mult
         * */
        void 
        mul_perm_mult (item_t mult);

        /**
         * \brief  Sets half of fracture length
         * \param  half_length Half of fracte length
         * */
        void 
        set_half_length (item_t half_length);

        /**
         * \brief  Sets angle between fracture and positive X direction
         * \param  theta An angle
         * */
        void 
        set_theta (item_t theta);

        /**
         * \brief  Sets skin factor
         * \param  skin
         * */
        void 
        set_skin (item_t skin);

        /**
         * \brief  Sets perforation (connection) status
         * \param  connection_status
         * */
        void 
        set_status (connection_status_type connection_status);

        /**
         * \brief  Sets transmissibility factor
         * \param  factor
         * */
        void 
        set_factor (item_t factor);

        /**
         * \brief  Sets wellbore diameter
         * \param  diameter
         * */
        void 
        set_diameter (item_t diameter);

        /**
         * \brief  Sets effective Kh 
         *         (permeability * thikness)
         * \param  Kh
         * */
        void 
        set_Kh (item_t kh);

        /**
         * \brief  Sets direction in which the well 
         *         penetrates the grid block
         * \param  
         * \return 
         * */
        void 
        set_direction (connection_direction_type direction);

        /**
         * \brief  Sets coordinates and index of grid block
         * \param  i
         * \param  j
         * \param  k
         * \param  n_block Number of grid block
         * \return 
         * */
        void 
        set_coord (index_t i, index_t j, index_t k, index_t n_block);

        /**
         * \brief  Sets depth
         * \param  mesh
         * */
        void 
        set_connection_depth (const sp_mesh_iface_t &mesh);

        /**
         * \brief  Returns index of grid block to
         *         which perforation belongs
         * \return Index of grid block
         * */
        index_t 
        n_block () const;

        /**
         * \brief  Returns connection_factor
         * \return Connection factor
         * */
        item_t 
        get_fact () const;

        /**
         * \brief  Sets bulkp
         * \param  bulkp
         * */
        void 
        set_bulkp (item_t bulkp);

        /**
         * \brief  Sets rate
         * \param  rate
         * \todo   Obsolete, should be removed 
         * */
        void 
        set_rate (item_t rate);

        /**
         * \brief  Sets head_term (?)
         * \param  head_term
         * */
        void 
        set_head_term (item_t head_term);

        /**
         * \brief  Sets current connection BHP
         * \param  cur_bhp
         * */
        void 
        set_cur_bhp (item_t cur_bhp);

        /**
         * \brief  Sets multiplier
         * \param  mult
         * */
        void 
        set_mult (item_t mult);

        /**
         * \brief  Sets number of segment which
         *         contains this perforation (connection)
         * \param  seg
         * */
        void 
        set_seg_number (const index_t &seg);

        /**
         * \brief  Is perforation (connection) is shut?
         * \return True if status of connection is shut
         * */
        bool 
        is_shut () const;

        /**
         * \brief  Returns status of connection
         * \return Connection status
         * */
        connection_status_type 
        get_status () const
          {
            return status_;
          }

        /**
         * \brief  Returns direction of connection
         * \return Connection direction
         * */
        connection_direction_type 
        get_dir () const
          {
            return dir_;
          }

        //! Returns current BHP
        item_t                
        get_cur_bhp () const;

        //! Returns density
        item_t
        get_density () const;

        //! Returns depth
        item_t
        get_connection_depth () const;

        //! Returns bulkp
        item_t
        get_bulkp () const;

        //! Returns i coord
        index_t
        i_coord () const;

        //! Return j coord
        index_t
        j_coord () const;

        //! Return k coord
        index_t
        k_coord () const;

        //! Return multiplier
        item_t
        mult () const;

        //! Returns head_term (?)
        item_t
        get_head_term () const;

        //! Returns segment number
        index_t
        get_seg_number () const;

        /**
         * \brief  Clears data
         * */
        virtual void 
        clear_data ();

        /**
         * \brief  Returns RW value
         * \todo   Obsolete, should be removed
         * */
        virtual array_ext <item_t> 
        get_rw_value   ();

        /**
         * \brief  Returns WR value
         * \todo   Obsolete, should be removed
         * */
        virtual array_ext <item_t> 
        get_wr_value   ();

        /**
         * \brief  Returns RR value
         * \todo   Obsolete, should be removed
         * */
        virtual array_ext <item_t> 
        get_rr_value   ();

        /**
         * \brief  Returns PS value
         * \todo   Obsolete, should be removed
         * */
        virtual array_ext <item_t> 
        get_ps_value   ();

        /**
         * \brief  Returns rate array
         * \todo   Obsolete, should be removed
         * */
        virtual array_ext <rhs_item_t> 
        get_rate_value ();

        /**
         * \brief  Returns rate data
         * \return Rate data
         * */
        const rate_data_t &
        rate () const
        {
          return rate_;
        }

        /**
         * \brief  Returns production part of rate data
         * \return Production part of rate data
         * */
        const rate_data_inner_t &
        rate_prod () const
        {
          return rate_.prod;
        }

        /**
         * \brief  Returns injection part of rate data
         * \return Injection part of rate data
         * */
        const rate_data_inner_t &
        rate_inj () const
        {
          return rate_.inj;
        }

      public:
        auto_value <item_t>           head_term;
        auto_value <item_t>           cur_bhp;
        auto_value <item_t>           connection_depth;
        auto_value <item_t>           density;
        auto_value <item_t>           bulkp;

        rate_data_t                   rate_;
        rate_data_t                   rate_rc_;

      private:

        template <typename strategy_t_x>
        friend struct wells::compute_factors::peaceman_model;

        template <typename strategy_t_x>
        friend struct wells::compute_factors::baby_odeh_model;

        auto_value <index_t, -1>      i_coord_;                     //!< i coord of connection
        auto_value <index_t, -1>      j_coord_;                     //!< j coord of connection
        auto_value <index_t, -1>      k_coord_;                     //!< k coord of connection

        auto_value <connection_status_type, connection_open>
        status_;                      //!< connection status

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
        //! blue-sky type declaration
        BLUE_SKY_TYPE_DECL_T (connection <strategy_t>);
      };


  } // namespace wells
} // namespace blue_sky


#endif  // #ifndef BS_WELL_CONNECTION_H_
