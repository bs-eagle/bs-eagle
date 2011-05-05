#ifndef INIT_MODEL_IFACE_HPP_04ee155e_754c_11e0_b98f_a3a6a578db4a
#define INIT_MODEL_IFACE_HPP_04ee155e_754c_11e0_b98f_a3a6a578db4a
/**
 *       \file  init_model_iface.hpp
 *      \brief  interface class for init_model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  03.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "bs_object_base.h"

namespace blue_sky 
{
  class idata;
  class rs_mesh_iface;
  class calc_model;
  class BS_API_PLUGIN init_model_iface : public objbase
  {
  public:

    virtual void
    init (BS_SP (calc_model) model, BS_SP (idata) data, BS_SP (rs_mesh_iface) mesh) = 0;

  };

}

#endif //

