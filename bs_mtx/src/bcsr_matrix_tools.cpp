/**
 * @file bcsr_matrix_tools.cpp
 * @brief
 * @date 2009-11-24
 */
#include "bcsr_matrix_tools.h"

#include <time.h>
#include <stdlib.h>
#include <map>

//#include "strategies.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  bcsr_matrix_tools::bcsr_matrix_tools (bs_type_ctor_param)
        : bcsr_matrix_tools_iface ()
    {
    }
  bcsr_matrix_tools::bcsr_matrix_tools (const bcsr_matrix_tools & /*src*/) : bs_refcounter ()
     {
     }

/*!
  \brief read matrix from ascii file
  \param matrix         -- reference to the BCSR matrix interface
  \param file_name -- name of the file
  \return 0 if success
*/
int
bcsr_matrix_tools::ascii_read_from_csr_format (sp_bcsr_t matrix, const std::string &file_name) const
{
  //locale_keeper lkeeper ("C", LC_NUMERIC);

  FILE *fp = 0;
  char buf[4096];
  int state = 0;
  t_long nc, nr, nnz = 0, nbs, b_sqr = 0;
  t_long row_ind = 0;
  t_long j_ind = 0, j = 0;
  t_long *rows_ptr;
  t_long *cols_ind;
  t_float *values;
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
          t_long n_rows = matrix->get_n_rows ();

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

/*!
  \brief write matrix to ascii file
  \param matrix         -- reference to the BCSR matrix interface
  \param file_name -- name of the file
  \return 0 if success
*/
int
bcsr_matrix_tools::ascii_write_to_csr_format (const sp_bcsr_t matrix,
                                              const std::string &file_name,
                                              const bool sort_cols) const
{
  FILE *fp = 0;
  t_long i, j, j1, j2, jj, jj1, jj2, counter, b_sqr = 0;
  t_long n_memory_sort_index = 10, n_row_cols, r_code = 0;
  spv_long sp_sort_index = BS_KERNEL.create_object (v_long::bs_type ());
  BS_ASSERT (sp_sort_index);
  t_long *sort_index = 0;

  t_long *rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
  t_long *cols_ind = &(*(matrix->get_cols_ind ()))[0];
  t_float *values = &(*(matrix->get_values ()))[0];
  t_long n_rows = matrix->get_n_rows ();
  t_long n_cols = matrix->get_n_cols ();
  t_long n_block_size = matrix->get_n_block_size ();
  t_long n_nnz = matrix->get_n_non_zeros ();
  b_sqr = n_block_size * n_block_size;

  fp = fopen (file_name.c_str (), "w");
  if (!fp)
    //TODO: write error message
    return -1;
  fprintf (fp, "// N_ROWS\tN_COLS\tN_NON_ZEROS\tN_BLOCK_SIZE\n");
  fprintf (fp, "%d\t%d\t%d\t%d\n", n_rows, n_cols, n_nnz, n_block_size);

  fprintf (fp, "// Rows indexes[1..n_rows] (with out 0)\n");
  for (i = 1; i <= n_rows; ++i)
    fprintf (fp, "%d\n", rows_ptr[i]);

  fprintf (fp, "// END of Rows indexes\n");


  fprintf (fp, "// Values n_non_zeros elements\n");
  fprintf (fp, "//COLUMN\tVALUE\n");
  for (i = 0, counter = 0; i < n_rows; ++i)
    {
      fprintf (fp, "// ROW %d\n", i);
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];

      if (sort_cols)
        {
          sp_sort_index->resize (j2 - j1);
          sort_index = &(*(sp_sort_index))[0];
          for (j = j1, jj = 0; j < j2; ++j, ++jj)
            {
              sort_index[jj] = j;
            }

          n_row_cols = j2 - j1;
          for (j = 0; j < n_row_cols - 1; ++j)
            {
              for (jj = j; jj < n_row_cols; ++jj)
                {
                  if (cols_ind[sort_index[j]] > cols_ind[sort_index[jj]])
                    {
                      jj1 = sort_index[j];
                      sort_index[j] = sort_index[jj];
                      sort_index[jj] = jj1;
                    }
                }
            }
          for (j = 0; j < n_row_cols; ++j, ++counter)
            {
              fprintf (fp, "%d", cols_ind[sort_index[j]]);
              jj1 = sort_index[j] * b_sqr;
              jj2 = jj1 + b_sqr;
                for (jj = jj1; jj < jj2; ++jj)
                  fprintf (fp, "\t%le", values[jj]);
              fprintf (fp, "\n");
            }
        }
      else
        {
          for (j = j1; j < j2; ++j, ++counter)
            {
              fprintf (fp, "%d", cols_ind[j]);
              jj1 = j * b_sqr;
              jj2 = jj1 + b_sqr;
                for (jj = jj1; jj < jj2; ++jj)
                  fprintf (fp, "\t%le", values[jj]);
              fprintf (fp, "\n");
            }
        }
    }
  if (counter != n_nnz)
    return -1;
  fprintf (fp, "// END OF FILE\n");
  fclose (fp);

  return 0;
}

int
bcsr_matrix_tools::random_init (sp_bcsr_t matrix,
                                const t_long new_n_rows,
                                const t_long new_n_block_size,
                                const t_double rand_value_dispersion,
                                const t_long elems_in_row
                               ) const
{
  int r_code = 0;
  t_long i, j, j1, j2, b_sqr;
  t_long n_non_zeros;


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

  t_float *values = &(*(matrix->get_values ()))[0];
  t_long *rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
  t_long *cols_ind = &(*(matrix->get_cols_ind ()))[0];

  srand ((unsigned)time( NULL ));

  for (i = 0; i < n_non_zeros * b_sqr; ++i)
    values[i] = (t_float)rand () / (t_float)RAND_MAX * (t_float)rand_value_dispersion;

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


int
bcsr_matrix_tools::gen_2d_laplas (sp_bcsr_t matrix, const t_long n) const
{
  int r_code = 0;
  t_long i, j;

  t_long elems_in_row = 5;
  t_long new_n_block_size = 1;

  // setup integer vars
  t_long N = n * n;
  t_long n_non_zeros = N * elems_in_row; // approx value, final nonzeros value can be equal or less

  r_code = matrix->init (N, N, new_n_block_size, n_non_zeros);

  if (r_code)
    //TODO: print error message
    return -2;

  t_float *values = &(*(matrix->get_values ()))[0];
  t_long *rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
  t_long *cols_ind = &(*(matrix->get_cols_ind ()))[0];

  rows_ptr[0] = 0;
  for (i = 0, j = 0; i < N; ++i)
    {
      values[j] = 4;
      cols_ind[j] = i;
      ++j;

      if (i - n >= 0)
        {
          values[j] = -1;
          cols_ind[j] = i - n;
          ++j;
        }

      if (i % n)
        {
          values[j] = -1;
          cols_ind[j] = i - 1;
          ++j;
        }

      if ((i + 1) % n)
        {
          values[j] = -1;
          cols_ind[j] = i + 1;
          ++j;
        }

      if (i + n < N)
        {
          values[j] = -1;
          cols_ind[j] = i + n;
          ++j;
        }

      rows_ptr[i + 1] = j;
    }

  //n_non_zeros = rows_ptr[new_n_rows];
  return 0;
}


int
bcsr_matrix_tools::gen_diagonal (sp_bcsr_t matrix, const t_long n,
                                 const t_long nb, const t_double val) const
{
  int r_code = 0;
  t_long i, k;
  t_long b_sqr = nb * nb;

  r_code = matrix->init (n, n, nb, n * b_sqr);
  if (r_code)
    //TODO: print error message
    return -2;

  t_float *values = &(*(matrix->get_values ()))[0];
  t_long *rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
  t_long *cols_ind = &(*(matrix->get_cols_ind ()))[0];

  rows_ptr[0] = 0;
  for (i = 0; i < n; ++i)
    {
      for (k = 0; k < b_sqr; ++k)
        {
          if (k % (nb + 1) == 0)
            values[i * b_sqr + k] = val;
          else
            values[i * b_sqr + k] = 0.0;
        }
      cols_ind[i] = i;
      rows_ptr[i + 1] = i + 1;
    }

  return 0;
}

int
bcsr_matrix_tools::dense_init (sp_bcsr_t matrix,
                               const t_long n_rows,
                               const t_long block_size,
                               const t_double rand_value_dispersion) const
{
  int r_code = 0;
  t_long b_sqr = block_size * block_size;

  srand ((unsigned)time( NULL ));

  r_code = matrix->init (n_rows, n_rows, block_size, n_rows * n_rows);
  if (r_code)
    return -2;

  t_float *values = &(*(matrix->get_values ()))[0];
  t_long *rows_ptr = &(*(matrix->get_rows_ptr ()))[0];
  t_long *cols_ind = &(*(matrix->get_cols_ind ()))[0];
  rows_ptr[0] = 0;


  for (t_long i = 0; i < n_rows * n_rows * b_sqr; ++i)
    {
      values[i] = (t_float)rand () / (t_float)RAND_MAX * (t_float)rand_value_dispersion;
    }
  for (t_long i = 0; i < n_rows; ++i)
    {
      t_long cl;
      t_long j1 = rows_ptr[i];

      rows_ptr[i + 1] = j1 + n_rows;

      // diagonal should be the first
      cols_ind[j1] = i;
      cl = j1 + 1;
      for (t_long j = 0; j < n_rows; ++j)
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

  BLUE_SKY_TYPE_STD_CREATE (bcsr_matrix_tools);
  BLUE_SKY_TYPE_STD_COPY (bcsr_matrix_tools);

  BLUE_SKY_TYPE_IMPL(bcsr_matrix_tools,  bcsr_matrix_tools_iface, "bcsr_matrix_tools", "Tools for Block CSR Matrix class", "Tools realization of Block CSR Matricies");
}  // blue_sky namespace
