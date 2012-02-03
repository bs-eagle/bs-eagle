/**
 *       \file  data_dimens.h
 *      \brief  data_dimens class declaration
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.11.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_EAGLE_DATA_DIMENS_H_
#define BS_EAGLE_DATA_DIMENS_H_

#include "conf.h"

namespace blue_sky {

  /**
   * \class data_dimens
   * \brief Holds X, Y, Z dimensions of data
   * */
  struct data_dimens
  {
    t_long nx;     //!< X dimension
    t_long ny;     //!< Y dimension
    t_long nz;     //!< Z dimension

    data_dimens (t_long nx = 0, t_long ny = 0, t_long nz = 0)
    : nx (nx)
    , ny (ny)
    , nz (nz)
    {
    }
  };


  namespace python {

    /**
     * \brief  Exports data_dimens to python
     * */
    void
    export_data_dimens ();

  } // namespace python

}  // namespace blue_sky


#endif // #ifndef BS_EAGLE_DATA_DIMENS_H_

