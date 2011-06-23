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

namespace blue_sky {

  /**
   * \class data_dimens
   * \brief Holds X, Y, Z dimensions of data
   * */
  struct data_dimens
  {
    long long nx;     //!< X dimension
    long long ny;     //!< Y dimension
    long long nz;     //!< Z dimension

    data_dimens (long long nx = 0, long long ny = 0, long long nz = 0)
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

