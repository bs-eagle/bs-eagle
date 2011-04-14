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
  class bcsr_matrix_tools_iface: public bs_node
    {

    public:
      typedef bcsr_matrix_iface                                 bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;

    public:
      //! destructor
      virtual ~bcsr_matrix_tools_iface ()
        {};

      // read matrix from ascii file
      virtual int ascii_read_from_csr_format (sp_bcsr_t matrix, const std::string &file_name) const = 0;

      virtual int random_init (sp_bcsr_t matrix,
                               const t_long new_n_rows,
                               const t_long new_n_block_size,
                               const t_double rand_value_dispersion,
                               const t_long elems_in_row
                               ) const = 0;
      virtual int gen_2d_laplas (sp_bcsr_t matrix, const t_long n) const = 0;
      virtual int dense_init (sp_bcsr_t matrix, const t_long n_rows, const t_long block_size,
                              const t_double rand_value_dispersion) const = 0;
      //virtual int gen_2d_laplas (bcsr_t &matrix, const int n) = 0;
    };


}//namespace blue_sky
#endif //__BCSR_MATRIX_TOOLS_IFACE_H


