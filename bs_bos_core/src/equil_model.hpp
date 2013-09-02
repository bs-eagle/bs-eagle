#ifndef EQUIL_MODEL_HPP_ffbe720e_755a_11e0_bdd6_67da5c425f32
#define EQUIL_MODEL_HPP_ffbe720e_755a_11e0_bdd6_67da5c425f32
/**
 *       \file  equil_model.hpp
 *      \brief  EQUIL Initialization model
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
  class BS_API_PLUGIN equil_model : public init_model_iface
  {
  public:

    BLUE_SKY_TYPE_DECL (equil_model);

    void
    init (BS_SP (calc_model) model, BS_SP (idata) data, BS_SP (rs_mesh_iface) mesh);
  };

}

#endif //

