/**
 *       \file  wells_compute_connection_factors.h
 *      \brief  Computes connection factors
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  07.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_WELLS_COMPUTE_CONNECTION_FACTORS_H_
#define BS_WELLS_COMPUTE_CONNECTION_FACTORS_H_

//#include "calc_well.h"

namespace blue_sky
  {
  namespace wells
    {
      class connection;
    namespace compute_factors
      {

        /**
         * \class peaceman_model
         * \brief Calculates connection factors (peaceman model)
         * */
      struct peaceman_model
        {

          typedef connection                      connection_t;
          typedef t_double                        item_t;
          typedef t_long                          index_t;
          typedef spv_double                      item_array_t;
          typedef smart_ptr <fi_params, true>     sp_params_t;

          typedef rs_mesh_iface                   mesh_iface_t;
          typedef smart_ptr <mesh_iface_t, true>  sp_mesh_iface_t;

          /**
           * \brief   computes connection factor by peaceman model.
           * \details d1 equals to msh->elements[k].d[dir_plane1]
           *          d2 equals to msh->elements[k].d[dir_plane2]
           *          d3 or d_ort equals to msh->elements[k].d[dir_orthogonal]
           * \param   connection
           * \param   internal_constants
           * \param   params
           * \param   mesh
           * \param   perm
           * \param   ntg
           * \param   ro_calc_flag
           * */
          static void 
          compute (connection_t &connection,
                   const physical_constants &internal_constants,
                   const sp_params_t &params,
                   const sp_mesh_iface_t &mesh,
                   const stdv_float &perm,
                   const stdv_float &ntg,
                   bool ro_calc_flag,
                   item_t *completion_thickness = 0);

          static item_t compute_grp_pi_mult (connection_t &connection);
        };
      
      struct completion_connection_factor 
        {
          typedef connection                      connection_t;
          typedef t_double                        item_t;
          typedef t_long                          index_t;
          typedef spv_double                      item_array_t;
          typedef smart_ptr <fi_params, true>     sp_params_t;

          typedef rs_mesh_iface                   mesh_iface_t;
          typedef smart_ptr <mesh_iface_t, true>  sp_mesh_iface_t;

          static void 
          compute (connection_t &connection,
                   const physical_constants &internal_constants,
                   const sp_params_t &params,
                   const sp_mesh_iface_t &mesh,
                   const stdv_float &perm,
                   const stdv_float &ntg);
        };
      //struct baby_odeh_model
      //{
      //  static void compute (well::connection &connection, physical_constants *internal_constants,
      //                       item_t d1, item_t d2, item_t d3,
      //                       item_t perm1, item_t perm2, item_t perm3,
      //                       const well::item_array_t &ntg,
      //                       bool ro_calc_flag);

      //private:

      //  static item_t F_table (item_t L, item_t d, item_t y_mid, int arg_type);
      //};

    } // namespace compute_factors
  } // namespace wells
} // namespace blue_sky

#endif  // #ifndef BS_WELLS_COMPUTE_CONNECTION_FACTORS_H_

