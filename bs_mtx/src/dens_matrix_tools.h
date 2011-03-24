#ifndef __DENS_MATRIX_TOOLS_H
#define __DENS_MATRIX_TOOLS_H
/** 
 * @file dens_matrix_tools.h
 * @brief 
 * @date 2009-11-24
 */
#include "bs_assert.h"
#include "bs_tree.h"

#include "dens_matrix_tools_iface.h"

#include <string>

namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  
  class dens_matrix_tools: public dens_matrix_tools_iface 
    {

    public:
      typedef dens_matrix_iface                                 dens_t;
      typedef smart_ptr<dens_t, true>                           sp_dens_t;

    public:
      BLUE_SKY_TYPE_DECL_T(dens_matrix_tools);
    public:
      //! destructor
      virtual ~dens_matrix_tools ()
        {};

      // read matrix from ascii file

      virtual int random_init (sp_dens_t matrix, 
                               const t_long ni, 
                               const t_long nj,
                               const t_long block_size,
                               const t_double rand_value_dispersion
                               ) const;
    };


}//namespace blue_sky
#endif //__DENS_MATRIX_TOOLS_IFACE_H

