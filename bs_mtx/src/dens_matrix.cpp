/** 
 * @file dens_matrix.cpp
 * @brief full matrix
 * @date 2009-12-08
 */
#include <sstream>

#include "bs_mtx_stdafx.h"
#include "dens_matrix.h" 
#include "pyublas/numpy.hpp"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  template <class strat_t>
  dens_matrix<strat_t>::dens_matrix (bs_type_ctor_param) 
        : dens_matrix_iface<strat_t> (),
        values (BS_KERNEL.create_object (fp_storage_array_t::bs_type ()))
    {
    }
  template <class strat_t>
   dens_matrix <strat_t>::dens_matrix (const this_t & /*src*/) : bs_refcounter (),
        values (BS_KERNEL.create_object (fp_storage_array_t::bs_type ()))
     {
     }

  // TODO: 
  template <class strat_t> int
  dens_matrix<strat_t>::matrix_vector_product (sp_fp_array_t v_, 
                                               sp_fp_array_t r_) const
    {
      i_type_t block_rows;
      i_type_t block_cols;
      i_type_t clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      if (n_cols != (i_type_t)v_->size () || n_rows != (i_type_t)r_->size ())
        return -1;

      fp_type_t *v = &(*v_)[0];
      fp_type_t *r = &(*r_)[0];

      block_rows = n_rows / clb;
      if (n_rows % clb)
        ++block_rows;
      
      block_cols = n_cols / clb;
      if (n_cols % clb)
        ++block_cols;

      for (i_type_t i = 0; i < block_rows; ++i)
        {
          for (i_type_t j = 0; j < block_cols; ++j)
            {
              block_mv_product (i, j, v, r);
            }
        }
      return 0;
    }

  // TODO: 
  template <class strat_t> int
  dens_matrix<strat_t>::matrix_vector_product_t (sp_fp_array_t v_, 
                                                 sp_fp_array_t r_) const
    {
      i_type_t block_rows;
      i_type_t block_cols;
      i_type_t clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      if (n_cols != (i_type_t)r_->size () || n_rows != (i_type_t)v_->size ())
        return -1;

      fp_type_t *v = &(*v_)[0];
      fp_type_t *r = &(*r_)[0];

      block_rows = n_rows / clb;
      if (n_rows % clb)
        ++block_rows;
      
      block_cols = n_cols / clb;
      if (n_cols % clb)
        ++block_cols;

      for (i_type_t i = 0; i < block_rows; ++i)
        {
          for (i_type_t j = 0; j < block_cols; ++j)
            {
              block_mv_product_t (i, j, v, r);
            }
        }
      return 0;
    }
          
  template <class strat_t> int
  dens_matrix<strat_t>::calc_lin_comb (fp_type_t alpha, 
                                       fp_type_t beta, 
                                       sp_fp_array_t u_, 
                                       sp_fp_array_t v_, 
                                       sp_fp_array_t r_) const
    {
      static const fp_type_t eps = fp_type_t (1.0e-12);
      i_type_t i;
      int r_code = 0;

      fp_type_t *v = &(*v_)[0];
      fp_type_t *r = &(*r_)[0];

      memset (r, 0, sizeof (fp_type_t) * n_rows);

      if (fabs (alpha) > eps)
        {
          r_code = matrix_vector_product (u_, r_);
          for (i = 0; i < n_rows; ++i)
            r[i] *= alpha;
        }
      if (fabs (beta) > eps)
        {
          for (i = 0; i < n_rows; ++i)
            r[i] += v[i] * beta;
        }
      return r_code;
    }

  template <class strat_t> int
  dens_matrix<strat_t>::init_by_matrix (sp_dens_matrix_iface_t matrix)
    {
      return init (matrix->get_n_rows (), matrix->get_n_cols (), matrix->get_calc_block_size ());
    }

  template <class strat_t> int
  dens_matrix<strat_t>::init (const i_type_t new_n_rows, 
                              const i_type_t new_n_cols, 
                              const i_type_t block_size)
    {
      npy_intp dims[2];
      if (new_n_rows < 1 || new_n_cols < 1)
        return -1;
      dims[0] = n_rows = new_n_rows;
      dims[1] = n_cols = new_n_cols;
      calc_block_size = block_size;
      values->resize (n_rows * n_cols);
      // TODO: add this method
      //values->numpy.reshape (2, dims);

      //std::fill (values->numpy.as_ublas ().data ().begin (), values->numpy.as_ublas ().data ().end (), 0);
      std::fill (values->begin (), values->end (), 0);
      return 0;
    }

  template <class strat_t> int
  dens_matrix<strat_t>::copy (sp_dens_matrix_iface_t matrix)
    {
      if (init (matrix->get_n_rows (), matrix->get_n_cols (), matrix->get_calc_block_size ()))
        return -3;
      values = matrix->get_values ()->clone ();
      //values->numpy.resize (matrix->get_values ()->numpy.size1 (), matrix->get_values ()->numpy.size2 ());
      //values->numpy.as_ublas () = matrix->get_values ()->numpy.as_ublas (); 
                              //                  matrix->get_values ()->numpy.as_ublas ().data ().end ());
       

      //memcpy (&values[0], &(matrix->get_values_const ())[0], n_rows * n_cols * sizeof (fp_storage_type_t));
      return 0;
    }

  template <class strat_t> void
  dens_matrix<strat_t>::block_mv_product (i_type_t row_block, 
                                          i_type_t col_block, 
                                          const fp_type_t *v, 
                                          fp_type_t *r) const
    {
      i_type_t block_n_rows;
      i_type_t block_n_cols;
      const fp_storage_type_t *block;
      const fp_type_t *block_v;
      fp_type_t *block_r;
      i_type_t clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      // calculate block start position
      //block = &values[row_block * clb * n_cols + col_block * clb];
      block = &(*values)[0] + row_block * clb * n_cols + col_block * clb;
      block_v = v + col_block * clb;
      block_r = r + row_block * clb;
      
      if ((row_block + 1) * clb <= n_rows)
        block_n_rows = clb;
      else
        block_n_rows = n_rows - row_block * clb;

      if ((col_block + 1) * clb <= n_cols)
        block_n_cols = clb;
      else
        block_n_cols = n_cols - col_block * clb;

      for (i_type_t i = 0; i < block_n_rows; ++i)
        {
          for (i_type_t j = 0; j < block_n_cols; ++j)
            {
              block_r[i] += block[i * n_rows + j] * block_v[j];
            }
        }
    }
  template <class strat_t> void
  dens_matrix<strat_t>::block_mv_product_t (i_type_t row_block, 
                                            i_type_t col_block, 
                                            const fp_type_t *v, 
                                            fp_type_t *r) const
    {
      i_type_t block_n_rows;
      i_type_t block_n_cols;
      const fp_storage_type_t *block;
      const fp_type_t *block_v;
      fp_type_t *block_r;
      i_type_t clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      // calculate block start position
      block = &(*values)[0] + row_block * clb * n_cols + col_block * clb;
      block_v = v + row_block * clb;
      block_r = r + col_block * clb;
      
      if ((row_block + 1) * clb <= n_rows)
        block_n_rows = clb;
      else
        block_n_rows = n_rows - row_block * clb;

      if ((col_block + 1) * clb <= n_cols)
        block_n_cols = clb;
      else
        block_n_cols = n_cols - col_block * clb;

      for (i_type_t i = 0; i < block_n_rows; ++i)
        {
          for (i_type_t j = 0; j < block_n_cols; ++j)
            {
              block_r[j] += block[i * n_rows + j] * block_v[i];
            }
        }
    }

#ifdef BSPY_EXPORTING_PLUGIN
  template <class strat_t> std::string
  dens_matrix<strat_t>::py_str () const
    {
      stringstream s;

      s << "Dens matrix\n--------------------------------\n";
      s << "\tRows:                  " << n_rows << "\n";
      s << "\tColumns:               " << n_cols << "\n";
      s << "\tCalc Block:            " << calc_block_size << "\n";
      s << "\tAllocated memory (Mb): " << get_allocated_memory_in_mbytes () << "\n";
      s << "--------------------------------\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(dens_matrix, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(dens_matrix, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix<base_strategy_fif>), 1, (dens_matrix_iface <base_strategy_fif> ),  "dens_matrix_fif",      "Dens Matrix class", "Realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix<base_strategy_did>), 1, (dens_matrix_iface <base_strategy_did> ),  "dens_matrix_did",      "Dens Matrix class", "Realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix<base_strategy_dif>), 1, (dens_matrix_iface <base_strategy_dif> ),  "dens_matrix_dif",      "Dens Matrix class", "Realization of Dens Matricies", false);

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix<base_strategy_flf>), 1, (dens_matrix_iface <base_strategy_flf> ),  "dens_matrix_flf",      "Dens Matrix class", "Realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix<base_strategy_dld>), 1, (dens_matrix_iface <base_strategy_dld> ),  "dens_matrix_dld",      "Dens Matrix class", "Realization of Dens Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (dens_matrix<base_strategy_dlf>), 1, (dens_matrix_iface <base_strategy_dlf> ),  "dens_matrix_dlf",      "Dens Matrix class", "Realization of Dens Matricies", false);
}  // blue_sky namespace
