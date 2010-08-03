#ifndef __DENS_MATRIX_TOOLS_H
#define __DENS_MATRIX_TOOLS_H
/** 
 * @file dens_matrix_tools.h
 * @brief 
 * @date 2009-11-24
 */
#include <string>

#include "bs_assert.h"
#include "bs_tree.h"

#include "dens_matrix_tools_iface.h"


namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  template <class strat_t>
  class dens_matrix_tools: public dens_matrix_tools_iface <strat_t>
    {

    public:
      typedef typename strat_t::fp_vector_type                  fp_vector_type_t;
      typedef typename strat_t::i_vector_type                   i_vector_type_t;
      typedef typename strat_t::fp_storage_vector_type          fp_storage_vector_type_t;
      typedef typename strat_t::fp_type_t                       fp_type_t;
      typedef typename strat_t::i_type_t                        i_type_t;
      typedef typename strat_t::fp_storage_type_t               fp_storage_type_t;
      typedef dens_matrix_iface <strat_t>                       dens_t;
      typedef smart_ptr<dens_t, true>                           sp_dens_t;

    public:
      BLUE_SKY_TYPE_DECL_T(dens_matrix_tools);
    public:
      //! destructor
      virtual ~dens_matrix_tools ()
        {};

      // read matrix from ascii file

      virtual int random_init (sp_dens_t matrix, 
                               const i_type_t ni, 
                               const i_type_t nj,
                               const i_type_t block_size,
                               const fp_type_t rand_value_dispersion
                               ) const;
    };


}//namespace blue_sky
#endif //__DENS_MATRIX_TOOLS_IFACE_H

