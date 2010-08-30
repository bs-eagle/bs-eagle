#ifndef __BCSR_MATRIX_TOOLS_IFACE_H
#define __BCSR_MATRIX_TOOLS_IFACE_H
/** 
 * @file bcsr_matrix_tools_iface.h
 * @brief 
 * @date 2009-11-24
 */

#include <string>

#include "bs_assert.h"
#include "bs_tree.h"

#include "bcsr_matrix_iface.h"


namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  template <class strat_t>
  class bcsr_matrix_tools_iface: public bs_node
    {

    public:
      //typedef typename strat_t::fp_vector_type                  fp_vector_type_t;
      //typedef typename strat_t::i_vector_type                   i_vector_type_t;
      //typedef typename strat_t::fp_storage_vector_type          fp_storage_vector_type_t;
      typedef typename strat_t::fp_type_t                       fp_type_t;
      typedef typename strat_t::i_type_t                        i_type_t;
      typedef typename strat_t::fp_storage_type_t               fp_storage_type_t;

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;

      typedef bcsr_matrix_iface<strat_t>                        bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;

    public:
      //! destructor
      virtual ~bcsr_matrix_tools_iface ()
        {};

      // read matrix from ascii file
      virtual int ascii_read_from_csr_format (sp_bcsr_t matrix, const std::string &file_name) const = 0;

      virtual int random_init (sp_bcsr_t matrix, 
                               const i_type_t new_n_rows, 
                               const i_type_t new_n_block_size,
                               const fp_type_t rand_value_dispersion, 
                               const i_type_t elems_in_row
                               ) const = 0;

      virtual int dense_init (sp_bcsr_t matrix, const i_type_t n_rows, const i_type_t block_size,
                              const fp_type_t rand_value_dispersion) const = 0;
      //virtual int gen_2d_laplas (bcsr_t &matrix, const int n) = 0;
    };


}//namespace blue_sky
#endif //__BCSR_MATRIX_TOOLS_IFACE_H


