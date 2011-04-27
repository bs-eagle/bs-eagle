/**
 *       \file  trans_multipliers.h
 *      \brief  Calculates transsmisibility mutlipliers
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  04.09.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete
 * */
#ifndef BS_TRANS_MULTIPLIERS_HPP_
#define BS_TRANS_MULTIPLIERS_HPP_

namespace blue_sky
  {

  struct trans_multipliers_calc
    {
      //trans_multipliers_calc (const sp_mesh_iface_t &mesh, const item_array_t &trans_mult, item_array_t &planes_trans)
      //: mesh (mesh)
      //, trans_mult (trans_mult)
      //, planes_trans (planes_trans)
      //{
      //}

      trans_multipliers_calc ()
      {
      }

      void apply ()
      {
        BS_ASSERT (false && "NOT IMPL YET");
      }
      void unapply ()
      {
        BS_ASSERT (false && "NOT IMPL YET");
      }

private:
      //const sp_mesh_iface_t     &mesh;
      //const item_array_t  &trans_mult;
      //item_array_t        &planes_trans;
    };

} // namespace blue_sky

#endif  // #ifndef BS_RANS_MULTIPLIERS_HPP_

