/** 
 * @file bcsr_matrix_tools.cpp
 * @brief 
 * @date 2009-11-24
 */
#include <time.h>
#include <stdlib.h>
#include <map>

#include "bs_matrix_stdafx.h"
#include "bcsr_matrix_tools.h"

#include "strategies.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  template <class strat_t>
  bcsr_matrix_tools<strat_t>::bcsr_matrix_tools (bs_type_ctor_param) 
        : bcsr_matrix_tools_iface<strat_t> ()
    {
    }
  template <class strat_t>
   bcsr_matrix_tools <strat_t>::bcsr_matrix_tools (const bcsr_matrix_tools <strat_t>& /*src*/) : bs_refcounter ()
     {
     }

/*!
  \brief read matrix from ascii file
  \param matrix         -- reference to the BCSR matrix interface
  \param file_name -- name of the file
  \return 0 if success
*/
template <class strat_t> int
bcsr_matrix_tools<strat_t>::ascii_read_from_csr_format (sp_bcsr_t matrix, const std::string &file_name) const 
{
  //locale_keeper lkeeper ("C", LC_NUMERIC);

  FILE *fp = 0;
  char buf[4096];
  int state = 0;
  i_type_t nc, nr, nnz = 0, nbs, b_sqr = 0;
  i_type_t row_ind = 0;
  i_type_t j_ind = 0, j = 0;
  i_type_t *rows_ptr;
  i_type_t *cols_ind;
  fp_storage_type_t *values;
  char *start_ptr, *end_ptr;

  fp = fopen (file_name.c_str (), "r");
  //BS_ASSERT (fp) (file_name);
  if (!fp)
    {
      //bs_throw_exception ("Can't open file");
      return -1;
    }

  while (fgets(buf, 4096, fp))
    {
      // skip comments
      if (!strncmp (buf, "//", 2))
        continue;
      if (state == 0) // read n_rows, n_cols, n_non_zeros
        {
          start_ptr = buf;
          nr = strtol (start_ptr, &end_ptr, 10);
          if (start_ptr == end_ptr)
            continue;
          start_ptr = end_ptr;
          nc = strtol (start_ptr, &end_ptr, 10);
          if (start_ptr == end_ptr)
            continue;
          start_ptr = end_ptr;
          nnz = strtol (start_ptr, &end_ptr, 10);
          if (start_ptr == end_ptr)
            continue;
          start_ptr = end_ptr;
          nbs = strtol (start_ptr, &end_ptr, 10);
          if (start_ptr == end_ptr)
            continue;

          if (matrix->init (nr, nc, nbs, nnz))
            {
              //bs_throw_exception ("Can't init matrix");
              return -45;
            }

          //BOSOUT << "nr = " << nr << ", nc = " << nc << ", nbs = " << nbs << ", nnz = " << nnz << bs_end;

          state = 1;

          rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
          cols_ind = &(*(matrix->get_cols_ind ()))[0];
          values = &(*(matrix->get_values ()))[0];

          rows_ptr[0] = 0;
          b_sqr = nbs * nbs;
        }
      else if (state == 1) // read rows_ptr
        {
          i_type_t n_rows = matrix->get_n_rows ();

          start_ptr = buf;
          rows_ptr[row_ind + 1] = strtol (start_ptr, &end_ptr, 10);
          if (start_ptr == end_ptr)
            continue;
          ++row_ind;
          if (row_ind == n_rows)
            state = 2;
        }
      else if (state == 2) // read values
        {
          start_ptr = buf;
          cols_ind[j_ind] = strtol (start_ptr, &end_ptr, 10);
          if (start_ptr == end_ptr)
            continue;
          for (j = 0; j < b_sqr; ++j)
            {
              start_ptr = end_ptr;
              values[j_ind * b_sqr + j] = strtod (start_ptr, &end_ptr);
              if (start_ptr == end_ptr)
                continue;
            }
          j_ind++;
          if (j_ind == nnz)
            break;
        }
    }

  //BOSOUT << "nnz = " << nnz << bs_end;

  fclose (fp);
  return 0;
}

template <class strat_t> int
bcsr_matrix_tools<strat_t>::random_init (sp_bcsr_t matrix, 
                                         const i_type_t new_n_rows, 
                                         const i_type_t new_n_block_size,
                                         const fp_type_t rand_value_dispersion, 
                                         const i_type_t elems_in_row
                                         ) const
{
  int r_code = 0;
  i_type_t i, j, j1, j2, b_sqr;
  i_type_t n_non_zeros;


  // check input data
  if (elems_in_row <= 0)
    // internal error
    return -1;

  // setup integer vars
  n_non_zeros = new_n_rows * elems_in_row; // approx value, final nonzeros value can be equal or less
  b_sqr = new_n_block_size * new_n_block_size;

  r_code = matrix->init (new_n_rows, new_n_rows, new_n_block_size, n_non_zeros);

  if (r_code)
    //TODO: print error message
    return -2;

  fp_storage_type_t *values = &(*(matrix->get_values ()))[0];
  i_type_t *rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
  i_type_t *cols_ind = &(*(matrix->get_cols_ind ()))[0];

  srand ((unsigned)time( NULL ));

  for (i = 0; i < n_non_zeros * b_sqr; ++i)
    values[i] = (fp_storage_type_t)rand () / (fp_storage_type_t)RAND_MAX * (fp_storage_type_t)rand_value_dispersion;

  rows_ptr[0] = 0;
  for (i = 0; i < new_n_rows; ++i)
    {
      rows_ptr[i + 1] = rows_ptr[i] +  elems_in_row;
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];
      cols_ind[j1++] = i; // diagonal should be the first
      
      for (j = j1; j < j2; ++j)
        {
          
          int cl = rand () % (int)((double)new_n_rows / (double)(j2 - j1)) + (j - j1) * new_n_rows / (j2 - j1);
          if (cl == i)
            {
              --j;
              continue;
            }
          cols_ind[j] = cl;
        }
    }

  //n_non_zeros = rows_ptr[new_n_rows];
  return 0;
}
      
template <class strat_t> int
bcsr_matrix_tools<strat_t>::dense_init (sp_bcsr_t matrix, 
                                        const i_type_t n_rows, 
                                        const i_type_t block_size,
                                        const fp_type_t rand_value_dispersion) const
{
  int r_code = 0;
  i_type_t b_sqr = block_size * block_size;

  srand ((unsigned)time( NULL ));

  r_code = matrix->init (n_rows, n_rows, block_size, n_rows * n_rows);
  if (r_code)
    return -2;

  fp_storage_type_t *values = &(*(matrix->get_values ()))[0];
  i_type_t *rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
  i_type_t *cols_ind = &(*(matrix->get_cols_ind ()))[0];
  rows_ptr[0] = 0;


  for (i_type_t i = 0; i < n_rows * n_rows * b_sqr; ++i)
    {
      values[i] = (fp_storage_type_t)rand () / (fp_storage_type_t)RAND_MAX * (fp_storage_type_t)rand_value_dispersion;
    }
  for (i_type_t i = 0; i < n_rows; ++i)
    {
      i_type_t cl;
      i_type_t j1 = rows_ptr[i];

      rows_ptr[i + 1] = j1 + n_rows;
      
      // diagonal should be the first
      cols_ind[j1] = i;
      cl = j1 + 1;
      for (i_type_t j = 0; j < n_rows; ++j) 
        {
          if (j != i)
            {
              cols_ind[cl] = j;
              ++cl;
            }
        }
    }
  return 0;
}
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(bcsr_matrix_tools, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(bcsr_matrix_tools, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_matrix_tools<base_strategy_fif>), 1,  (bcsr_matrix_tools_iface <base_strategy_fif> ), "bcsr_matrix_tools_fif", "Tools for Block CSR Matrix class", "Tools realization of Block CSR Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_matrix_tools<base_strategy_did>), 1,  (bcsr_matrix_tools_iface <base_strategy_did> ), "bcsr_matrix_tools_did", "Tools for Block CSR Matrix class", "Tools realization of Block CSR Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_matrix_tools<base_strategy_dif>), 1,  (bcsr_matrix_tools_iface <base_strategy_dif> ), "bcsr_matrix_tools_dif", "Tools for Block CSR Matrix class", "Tools realization of Block CSR Matricies", false);

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_matrix_tools<base_strategy_flf>), 1,  (bcsr_matrix_tools_iface <base_strategy_flf> ), "bcsr_matrix_tools_flf", "Tools for Block CSR Matrix class", "Tools realization of Block CSR Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_matrix_tools<base_strategy_dld>), 1,  (bcsr_matrix_tools_iface <base_strategy_dld> ), "bcsr_matrix_tools_dld", "Tools for Block CSR Matrix class", "Tools realization of Block CSR Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_matrix_tools<base_strategy_dlf>), 1,  (bcsr_matrix_tools_iface <base_strategy_dlf> ), "bcsr_matrix_tools_dlf", "Tools for Block CSR Matrix class", "Tools realization of Block CSR Matricies", false);
}  // blue_sky namespace
