#ifndef __DENS_MATRIX_TOOLS_IFACE_H
#define __DENS_MATRIX_TOOLS_IFACE_H
/** 
 * @file dens_matrix_tools_iface.h
 * @brief 
 * @date 2009-11-24
 */

#include <string>

#include "bs_assert.h"
#include "bs_tree.h"

#include "dens_matrix_iface.h"


namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN dens_matrix_tools_iface: public bs_node
    {
    public:
      typedef dens_matrix_iface                                 dens_t;
      typedef smart_ptr<dens_t, true>                           sp_dens_t;

    public:
      //! destructor
      virtual ~dens_matrix_tools_iface ()
        {};

      virtual int random_init (sp_dens_t matrix, 
                               const t_long ni, 
                               const t_long nj,
                               const t_long block_size,
                               const t_double rand_value_dispersion 
                               ) const = 0;
    };
}//namespace blue_sky
#endif //__DENS_MATRIX_TOOLS_IFACE_H


