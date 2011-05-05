#ifndef EXPLICIT_MODEL_HPP_b3428474_754b_11e0_a55e_8fd394202658
#define EXPLICIT_MODEL_HPP_b3428474_754b_11e0_a55e_8fd394202658
/**
 *       \file  explicit_model.hpp
 *      \brief  EXPLICIT Initialization model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  03.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "init_model_iface.hpp"

namespace blue_sky 
{
  class idata;
  class rs_mesh_iface;
  class calc_model;
  class BS_API_PLUGIN explicit_model : public init_model_iface
  {
  public:

    BLUE_SKY_TYPE_DECL (explicit_model);

    void
    init (BS_SP (calc_model) model, BS_SP (idata) data, BS_SP (rs_mesh_iface) mesh);
  };

}

#endif //

