/*!
* \file csr_matrix.cpp
* \brief implementation of compressed sparse row matrix
* \author Borschuk Oleg
* \date 2006-06-16
*/
#ifdef _MPI
#include <mpi.h>
#endif //_MPI

#include <memory.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "csr_matrix.h"
//#include "debug_tools.h"
#include "memory_macroses.h"
//#include "debug_macroses.h"
#include "matrix_macroses.h"
//#include "omp_tools.h"
#include "pure_mesh.h"

//! default constructor
csr_matrix::csr_matrix ()
{
  n_memory_cols_ind = 0;
  n_memory_values = 0;
  n_memory_n_rows = 0;
  // storage
  values = 0;

  // matrix indexes
  rows_ptr = 0;
  cols_ind = 0;
  diag_ind = 0;
  type = MATRIX_TYPE_CSR;
  own_memory = 1;
}


//! default destructor
csr_matrix::~csr_matrix ()
{
  this->mat_base::~mat_base ();

  // dim vars
  n_memory_cols_ind = 0;
  n_memory_values = 0;
  n_memory_n_rows = 0;

  if (own_memory)
    {
      // storage
      FI_FREE_ARRAY (values);

      // matrix indexes
      FI_FREE_ARRAY (rows_ptr);
      FI_FREE_ARRAY (cols_ind);
    }

  FI_FREE_ARRAY (diag_ind);
}

//! init by parent matrix
int csr_matrix::init (csr_matrix *matrix)
{
  init (matrix->n_rows, matrix->n_cols, matrix->n_block_size, matrix->get_n_non_zeros ());
  memcpy (rows_ptr, matrix->get_rows_ptr (), (n_rows + 1) * sizeof (int));
  memcpy (cols_ind, matrix->get_cols_ind (), (get_n_non_zeros()) * sizeof (int));
  memcpy (values, matrix->get_values(), (get_n_non_zeros()) * n_block_size * n_block_size * sizeof (double));
  return 0;
}

void 
csr_matrix::setup_all (const int new_n_rows, const int new_n_block_size, double *new_values, int *new_rows_ptr, int *new_cols_ind)
{
  if (own_memory)
    {
      // storage
      FI_FREE_ARRAY (values);

      // matrix indexes
      FI_FREE_ARRAY (rows_ptr);
      FI_FREE_ARRAY (cols_ind);
      n_memory_n_rows = 0;
      n_memory_cols_ind = 0;
      n_memory_values = 0;
    }
  n_rows = new_n_rows;
  n_block_size = new_n_block_size;
  n_cols = n_rows;
  values = new_values;
  rows_ptr = new_rows_ptr;
  cols_ind = new_cols_ind;
  own_memory = 0;
}
#if 0
void
csr_matrix::set_data (double *new_values, int *new_rows_ptr, int *new_cols_ind)
{
  values = new_values;
  rows_ptr = new_rows_ptr;
  cols_ind = new_cols_ind;
  own_memory = 0;
}
#endif //0

/*!
  \brief initialize matrix
  \param new_n_rows -- number of matrix rows (in blocks)
  \param new_n_cols -- number of matrix columns (in blocks)
  \param new_n_blok_size -- size of block
  \param new_n_non_zeros -- number of non zeros element in matrix (in blocks)
  \return 0 if success
*/
int
csr_matrix::init (const int new_n_rows, const int new_n_cols, const int new_n_block_size,
                  const int new_n_non_zeros)
{
  if (new_n_block_size < 1)
    return -1;

  if (new_n_non_zeros < 0)
    return -2;

  init_struct (new_n_rows, new_n_cols, new_n_non_zeros);

  n_block_size = new_n_block_size;
  if (alloc_values (new_n_non_zeros))
    return -3;

  return 0;
}

/*!
  \brief initialize matrix struct
  \param new_n_rows -- number of matrix rows (in blocks)
  \param new_n_cols -- number of matrix columns (in blocks)
  \param new_n_non_zeros -- number of non zeros element in matrix (in blocks)
  \return 0 if success
*/
int
csr_matrix::init_struct (const int new_n_rows, const int new_n_cols, const int new_n_non_zeros)
{
  if (new_n_cols < 1)
    return -1;

  n_cols = new_n_cols;
  if (alloc_rows_ptr (new_n_rows) || alloc_cols_ind (new_n_non_zeros))
    return -2;
  memset (rows_ptr, 0, sizeof (int) * (new_n_rows + 1));
  if (n_cols == n_rows)
    is_square = 1;
  return 0;
}

int
csr_matrix::alloc_rows_ptr (const int new_n_rows)
{
  int r_code = 0;

  if (new_n_rows < 1)
    return -1;
  n_rows = new_n_rows;
  if (n_rows < n_memory_n_rows)
    return 0;
  n_memory_n_rows = n_rows;
  if (!own_memory)
    {
      rows_ptr = 0;
      cols_ind = 0;
      values = 0;
      own_memory = 1;
    }
  FI_INT_ARRAY_REALLOCATOR (rows_ptr, n_memory_n_rows + 1, r_code);
  //FI_INT_ARRAY_REALLOCATOR (diag_ind, n_memory_n_rows, r_code);

  if (r_code)
    //TODO: print error message
    return -2;
  return 0;
}

int
csr_matrix::alloc_cols_ind_and_values (const int new_n_non_zeros)
{
  if (!own_memory)
    return -1;
  if (alloc_cols_ind (new_n_non_zeros) || alloc_values (new_n_non_zeros))
    return -2;
  return 0;
}

int
csr_matrix::alloc_values (const int new_n_non_zeros)
{
  int r_code = 0;
  int b_sqr = n_block_size * n_block_size;
  int nnz = new_n_non_zeros;
  if (!own_memory)
    return -1;
  if (nnz < 1)
    nnz = 1;

  if (b_sqr * nnz > n_memory_values)
    {
      FI_DOUBLE_ARRAY_REALLOCATOR (values, b_sqr * nnz, r_code);
      memset (values, 0, b_sqr * nnz * sizeof (double));
      n_memory_values = b_sqr * nnz;
    }

  if (r_code)
    //TODO: print error message
    return -2;
  return 0;
}

int
csr_matrix::alloc_cols_ind (const int new_n_non_zeros)
{
  int r_code = 0;
  int nnz = new_n_non_zeros;
  if (!own_memory)
    return -1;
  if (nnz < 1)
    nnz = 1;

  if (nnz > n_memory_cols_ind)
    {
      FI_INT_ARRAY_REALLOCATOR (cols_ind, nnz, r_code);
      memset (cols_ind, -1, nnz * sizeof (int));
      n_memory_cols_ind = nnz;
    }
  if (r_code)
    //TODO: print error message
    return -2;
  return 0;
}


int
csr_matrix::rand_init (const int new_n_rows, const int new_n_cols,
                   const double rand_value_dispersion, int elems_in_row)
{
  int r_code = 0;
  int i, j, j1, j2;
  int n_non_zeros;

  // check input data
  if (elems_in_row <= 0)
    // internal error
    return -1;

  // setup integer vars
  n_non_zeros = new_n_rows * elems_in_row; // approx value, final nonzeros value can be equal or less

  r_code = init (new_n_rows, new_n_cols, 1, n_non_zeros);

  if (r_code)
    //TODO: print error message
    return -2;

  srand ((unsigned)time( NULL ));

  for (i = 0; i < n_non_zeros; i ++)
    values[i] = (double)rand () / (double)RAND_MAX * rand_value_dispersion;

  rows_ptr[0] = 0;
  for (i = 0; i < n_rows; i ++)
    {
      rows_ptr[i + 1] = rows_ptr[i] + rand () % (elems_in_row + 1);
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];
      for (j = j1; j < j2; j++)
      {
        cols_ind[j] = rand () % (int)((double)n_cols / (double)(j2-j1)) + (j - j1) * n_cols / (j2-j1);
      }
    }

  n_non_zeros = rows_ptr[n_rows];

  return 0;
}

/*!
  \brief matrix vector product
  \param v -- vector
  \param r -- result vector
  \return 0 if success
*/
int
csr_matrix::matrix_vector_product (double *v, double *r)
{
  int i, j;
  int j1, j2;
  int cl;
  double *v_block, *r_block, *m_block;
  int b_sqr = n_block_size * n_block_size;

 // OMP_TIME_MEASURE_START (csr_matrix_vector_product_timer);

  CH_PTR (values);
  CH_PTR (cols_ind);
  CH_PTR (rows_ptr);

  CH_PTR (v);
  CH_PTR (r);

  ////////////////////////////////////////////
  // loop through rows
  ////////////////////////////////////////////

#ifdef CSR_MATRIX_VECTOR_PRODUCT_PARALLEL
  #pragma omp parallel for private (r_block, j1, j2, j, cl, m_block, v_block)
#endif //CSR_MATRIX_VECTOR_PRODUCT_PARALLEL

  for (i = 0; i < n_rows; ++i)
    {
      r_block = r + i * n_block_size;
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];

      ////////////////////////////////////////////
      // loop through columns in row
      ////////////////////////////////////////////
      for (j = j1; j < j2; ++j)
        {
          cl = cols_ind[j];
          m_block = values + j * b_sqr;
          v_block = v + cl * n_block_size;
          MV_PROD (n_block_size, m_block, v_block, r_block);
        }
    }

//  OMP_TIME_MEASURE_END (csr_matrix_vector_product_timer);
  return 0;
}

int
csr_matrix::build_transpose (csr_matrix *matrix, int rows_offset, int cols_offset, int new_n_rows)
{
  int i,j;
  int row_ind, j1, j2;
  int nnz;

  nnz = matrix->get_n_non_zeros ();

  if (new_n_rows == 0)
    {
      new_n_rows = matrix->n_cols;
    }

  init (new_n_rows, matrix->n_rows, matrix->n_block_size, nnz);

  // Count the number of entries in each column of matrix (row of transposed matrix) and fill the rows_ptr array.

  for (i = 0; i < nnz; i++)
    {
      ++rows_ptr[matrix->cols_ind[i] + 1 - rows_offset];
    }

  for (i = 2; i <= n_rows; i++)
    {
      rows_ptr[i] += rows_ptr[i - 1];
    }

  // Load the values and column numbers of transposed matrix


  for (i = 0; i < matrix->n_rows; i++)
    {
      j1 = matrix->rows_ptr[i];
      j2 = matrix->rows_ptr[i + 1];
      for (j = j1; j < j2; j++)
        {
          row_ind = matrix->cols_ind[j] - rows_offset;
          cols_ind[rows_ptr[row_ind]] = i + cols_offset;
          values[rows_ptr[row_ind]] = matrix->values[j];
          rows_ptr[row_ind]++;
        }
      }
  // rows_ptr now points to the *end* of the jth row of entries
  // instead of the beginning.  Restore rows_ptr to front of row.

  for (i = n_rows; i > 0; i--)
    {
      rows_ptr[i] = rows_ptr[i-1];
    }

  if (rows_ptr)
    rows_ptr[0] = 0;

  return 0;
}

int
csr_matrix::build_transpose_struct (csr_matrix *matrix, int rows_offset, int cols_offset, int new_n_rows)
{
  int i,j;
  int row_ind, j1, j2;
  int nnz;

  nnz = matrix->get_n_non_zeros ();

  if (new_n_rows == 0)
    {
      new_n_rows = matrix->n_cols;
    }

  init (new_n_rows, matrix->n_rows, matrix->n_block_size, nnz);

  // Count the number of entries in each column of matrix (row of transposed matrix) and fill the rows_ptr array.

  for (i = 0; i < nnz; i++)
    {
      ++rows_ptr[matrix->cols_ind[i] + 1 - rows_offset];
    }

  for (i = 2; i <= n_rows; i++)
    {
      rows_ptr[i] += rows_ptr[i - 1];
    }

  // Load the values and column numbers of transposed matrix


  for (i = 0; i < matrix->n_rows; i++)
    {
      j1 = matrix->rows_ptr[i];
      j2 = matrix->rows_ptr[i + 1];
      for (j = j1; j < j2; j++)
        {
          row_ind = matrix->cols_ind[j] - rows_offset;
          cols_ind[rows_ptr[row_ind]] = i + cols_offset;
          rows_ptr[row_ind]++;
        }
      }
  // rows_ptr now points to the *end* of the jth row of entries
  // instead of the beginning.  Restore rows_ptr to front of row.

  for (i = n_rows; i > 0; i--)
    {
      rows_ptr[i] = rows_ptr[i-1];
    }

  rows_ptr[0] = 0;

  return 0;
}


/*!
 * \brief write matrix to file with given name
 *
 * \param file_name -- file to write
 *
 * \return 0 if success
 */
int
csr_matrix::write_matrix_to_file (const char *file_name, int sort_cols)
{
  FILE *fp = 0;
  int i, j, j1, j2, jj, jj1, jj2, counter, b_sqr = 0;
  int *sort_index = 0, n_memory_sort_index = 10, n_row_cols, r_code = 0;

  CH_PTR (file_name);
  CH_PTR (cols_ind);
  CH_PTR (rows_ptr);
  CH_PTR (values);

  b_sqr = n_block_size * n_block_size;
  if (sort_cols)
    {
      FI_INT_ARRAY_ALLOCATOR (sort_index, n_memory_sort_index, r_code);
      if (r_code)
        //TODO: write error message
        return -1;
    }
  fp = fopen (file_name, "w");
  if (!fp)
    //TODO: write error message
    return -1;
  fprintf (fp, "// N_ROWS\tN_COLS\tN_NON_ZEROS\tN_BLOCK_SIZE\n");
  fprintf (fp, "%d\t%d\t%d\t%d\n", n_rows, n_cols, get_n_non_zeros (), n_block_size);

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

          while (j2 - j1 > n_memory_sort_index)
            {
              n_memory_sort_index *= 2;
              FI_INT_ARRAY_REALLOCATOR (sort_index, n_memory_sort_index, r_code);
              if (r_code)
                //TODO: write error message
                return -1;
            }

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
  if (counter != get_n_non_zeros ())
    return -1;
  fprintf (fp, "// END OF FILE\n");
  fclose (fp);
  FI_FREE_ARRAY (sort_index);
  return 0;
}
#ifdef _MPI
int
csr_matrix::mpi_write_matrix_to_file (const char *file_name, int sort_cols)
{
  char mpi_file_name[100];
  int rank, size;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  sprintf (mpi_file_name, "%s_%d_of_%d", file_name, rank, size);
  size = write_matrix_to_file (mpi_file_name, sort_cols);

  //delete mpi_file_name;
  return size;
}
#endif //_MPI


/** 
 * @brief read matrix from ascii format
 *        FORMAT:
 *              n_rows
 *              n_cols
 *              nnz
 *              rows[0]
 *              ...
 *              rows[n_rows]
 *              cols[0]
 *              ...
 *              cols[nnz - 1]
 *              values[0]
 *              ...
 *              values[nnz - 1]
 * 
 * @param file_name -- file name <INPUT>
 * 
 * @return 0 if success
 */
int 
csr_matrix::read_matrix_in_ascii_format (char *file_name)
{
  FILE *fp = 0;
  char buf[4096];
  int nc, nr, nnz = 0;
  int i;
  char *start_ptr, *end_ptr;
  

  CH_PTR (file_name);

  fp = fopen (file_name, "r");
  if (!fp)
    //TODO: write error message
    return -1;

  if (!fgets(buf, 4096, fp))
    return -1;
  start_ptr = buf;
  nr = strtol (start_ptr, &end_ptr, 10);

  if (!fgets(buf, 4096, fp))
    return -1;
  start_ptr = buf;
  nc = strtol (start_ptr, &end_ptr, 10);

  if (!fgets(buf, 4096, fp))
    return -1;
  start_ptr = buf;
  nnz = strtol (start_ptr, &end_ptr, 10);

  if (init (nr, nc, 1, nnz))
    return -45;
  printf ("READ %d %d %d\n", nr, nc, nnz);
  for (i = 0; i <= nr; ++i)
    {
      if (!fgets(buf, 4096, fp))
        return -1;
      start_ptr = buf;
      rows_ptr[i] = strtol (start_ptr, &end_ptr, 10);
    }
  printf ("ROWS FINiSHED\n");
  for (i = 0; i < nnz; ++i)
    {
      if (!fgets(buf, 4096, fp))
        return -1;
      start_ptr = buf;
      cols_ind[i] = strtol (start_ptr, &end_ptr, 10);
    }
  printf ("COLS FINiSHED\n");
  for (i = 0; i < nnz; ++i)
    {
      if (!fgets(buf, 4096, fp))
        return -1;
      start_ptr = buf;
      values[i] = strtod (start_ptr, &end_ptr);
    }
  printf ("VALUES FINiSHED\n");
  fclose (fp);
  return 0;
}
/*!
 * \brief read matrix from file with given name
 *
 * \param file_name -- file to read from
 *
 * \return 0 if suceess
 */
int
csr_matrix::read_matrix_from_file (const char *file_name)
{
  FILE *fp = 0;
  char buf[4096];
  int state = 0;
  int nc, nr, nnz = 0, nbs, b_sqr = 0;
  int row_ind = 0;
  int j_ind = 0, j = 0;
  char *start_ptr, *end_ptr;

  CH_PTR (file_name);

  fp = fopen (file_name, "r");
  if (!fp)
    //TODO: write error message
    return -1;

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

          if (init (nr, nc, nbs, nnz))
            return -45;
          state = 1;
          rows_ptr[0] = 0;
          b_sqr = nbs * nbs;
        }
      else if (state == 1) // read rows_ptr
        {
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
  fclose (fp);
  return 0;
}

int
csr_matrix::rand_init_symm (const int nx, const int ny, const int nz, const int nb,
                            const double rand_value_dispersion)
{

  int r_code = 0;
  int i, j, ind[6], k;
  int n, nnz, b_sqr, j_block;
  double d, *v_block;

  // check input data
  if (nx < 1 || ny < 1 || nz < 1)
    // internal error
    return -1;

  n = nx * ny * nz;
  nnz = n * 7 * nb; // approx value, final nonzeros value can be equal or less

  r_code = init (n, n, nb, nnz);
  b_sqr = n_block_size * n_block_size;

  if (r_code)
    //TODO: print error message
    return -2;

  srand ((unsigned)time( NULL ));

#ifdef _MPI_SEQ_RUN
  srand (1);
#endif

  for (i = 0; i < nnz; ++i)
    values[i] = (double)rand () / (double)RAND_MAX * rand_value_dispersion;

  int diag_j = 0;
  rows_ptr[0] = 0;
  for (i = 0, j = 0; i < n_rows; ++i)
    {
      ind[0] = i - nx * ny;
      ind[1] = i - nx;
      ind[2] = i - 1;
      ind[3] = i + 1;
      ind[4] = i + nx;
      ind[5] = i + nx * ny;

      d = 0;
      cols_ind[j] = i;// diagonal
      diag_j = j; // store diagonal position
      j++;

      for (k = 0; k < 6; ++k)// fill offdiagonal and add to sum
        if (ind[k] >= 0 && ind[k] < n)
          {
            cols_ind[j] = ind[k];
            v_block = values + j * b_sqr;
            for (j_block = 0; j_block < b_sqr; j_block++)
              {
                v_block[j_block] = - (double)rand () / (double)RAND_MAX * rand_value_dispersion;
                d += v_block[j_block];
              }
            ++j;
          }

      // at last, fill diagonal values
      v_block = values + diag_j * b_sqr;
      for (j_block = 0; j_block < b_sqr; j_block++)
        {
          v_block[j_block] = -d + (double)rand () / (double)RAND_MAX * rand_value_dispersion;
        }

      rows_ptr[i + 1] = j;
    }
  return 0;
}

int
csr_matrix::gen_2d_laplas (const int n)
{
  int r_code = 0;
  int i, j;
  int N, nnz;

  // check input data
  if (n < 1)
    // internal error
    return -1;
  N = n * n;

  nnz = N * 5; // approx value, final nonzeros value can be equal or less

  r_code = init (N, N, 1, nnz);

  if (r_code)
    //TODO: print error message
    return -2;

  rows_ptr[0] = 0;
  for (i = 0, j = 0; i < n_rows; ++i)
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
  return 0;
}

int
csr_matrix::alloc_diag_ind ()
{
  int r_code = 0;
  FI_INT_ARRAY_REALLOCATOR (diag_ind, n_rows, r_code);
  if (r_code)
    return -2;
  return 0;
}

/*!
 * \brief initialize diagonal indexes
 *
 * \return 0 if success
 */
int
csr_matrix::init_diag_ind (int diag_may_not_exist, int row_offset)
{
  int i, j1, j2, j;
  int r_code = 0;

  CH_PTR (rows_ptr);
  CH_PTR (cols_ind);
  CH_PTR (values);
  if (n_rows != n_cols)
    return -1;
  if (n_rows < 1)
    return -1;

  FI_INT_ARRAY_REALLOCATOR (diag_ind, n_rows, r_code);
  if (r_code)
    return -2;

  for (i = 0; i < n_rows; ++i)
    {
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];
      for (j = j1; j < j2; ++j)
        {
          if (cols_ind[j] == i + row_offset)
            {
              diag_ind[i] = j;
              break;
            }
        }
      if (j == j2) // diagonal not found
        {
          diag_ind[i] = -1;
          if (!diag_may_not_exist)
            return -4;
        }
    }
  return 0;
}

/*!
 * \brief initialize next to diagonal indexes
 *
 * \return 0 if success
 */
int
csr_matrix::init_next_to_diag_ind (int row_offset, int *col_map)
{
  int i, j1, j2, j;
  int r_code = 0;

  CH_PTR (rows_ptr);
  CH_PTR (cols_ind);
  CH_PTR (values);
  if (n_rows < 1)
    return -1;

  FI_INT_ARRAY_REALLOCATOR (diag_ind, n_rows, r_code);
  if (r_code)
    return -2;

  if (col_map)
    {
      for (i = 0; i < n_rows; ++i)
        {
          j1 = rows_ptr[i];
          j2 = rows_ptr[i + 1];
          diag_ind[i] = j2;
          for (j = j1; j < j2; ++j)
            {
              if (col_map[cols_ind[j]] > (i + row_offset))
                {
                  diag_ind[i] = j;
                  break;
                }
            }
        }
    }
  else
    {
      for (i = 0; i < n_rows; ++i)
        {
          j1 = rows_ptr[i];
          j2 = rows_ptr[i + 1];
          diag_ind[i] = j2;
          for (j = j1; j < j2; ++j)
            {
              if (cols_ind[j] > (i + row_offset))
                {
                  diag_ind[i] = j;
                  break;
                }
            }
        }
    }
  return 0;
}


/*!
 * \brief print staticstic of matrix
 * \return 0 if success
 */
int
csr_matrix::print_statistic ()
{
  CH_PTR (values);
  CH_PTR (cols_ind);
  CH_PTR (rows_ptr);

  printf ("---------Matrix Statistic---------\n");
  printf ("-- NNZ       %d\n", rows_ptr[n_rows]);
  printf ("-- STENS     %d\n", rows_ptr[n_rows] / n_rows);
  printf ("----------------------------------\n");
  return 0;
}

//! return number of nonzeros elements
int
csr_matrix::get_n_non_zeros ()
{
  if (rows_ptr)
    return rows_ptr[n_rows];
  else
    return 0;
}

/*!
 * \brief calculate r += A^T * v
 *
 * \param v -- vector
 * \param r -- result vector
 *
 * \return 0 -- if success
 */
int
csr_matrix::matrix_vector_product_t (double *v, double *r)
{
  int i, j1, j2, j, cl;

  //OMP_TIME_MEASURE_START (csr_matrix_vector_product_t_timer);

  CH_PTR (v);
  CH_PTR (r);
  CH_PTR (values);
  CH_PTR (rows_ptr);
  CH_PTR (cols_ind);
  if (n_block_size > 1)
    return -1;

  ////////////////////////////////////////////
  // loop through columns in transform matrix
  ////////////////////////////////////////////
  for (i = 0; i < n_rows; ++i)
    {
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];

      ////////////////////////////////////////////
      // loop through rows in transform matrix
      ////////////////////////////////////////////
      for (j = j1; j < j2; ++j)
        {
          cl = cols_ind[j];
          r[cl] += v[i] * values[j];
        }
    }
  //OMP_TIME_MEASURE_END (csr_matrix_vector_product_t_timer);
  return 0;
}

/*!
 * \brief calculate r = alpha * A * u + beta * v
 *
 * \param alpha -- coef
 * \param beta  -- coef
 * \param u     -- vector
 * \param v     -- vector
 * \param r     -- result vector
 *
 * \return 0 if success
 */
int
csr_matrix::calc_lin_comb (const double alpha, const double beta,
                           double *u, double *v, double *r)
{
  const double eps = 1.0e-12;
  int i, n, r_code = 0;
  double d;
//  OMP_TIME_MEASURE_START (csr_calc_lin_comb_timer);
  CH_PTR (r);

  d = beta;
  if (fabs (alpha) > eps)
    {
      CH_PTR (u);
      CH_PTR (rows_ptr);
      CH_PTR (cols_ind);
      CH_PTR (values);
      d /= alpha;
    }
  if (fabs (beta) > eps)
    {
      CH_PTR (v);
    }
  n = n_rows * n_block_size;


  if (fabs (beta) > eps)
#ifdef CSR_CALC_LIN_COMB_PARALLEL
  #pragma omp parallel for
#endif //CSR_CALC_LIN_COMB_PARALLEL
    for (i = 0; i < n; ++i)
      r[i] = v[i] * d;
  else
    memset (r, 0, sizeof (double) * n);

  if (fabs (alpha) > eps)
    {
      r_code = matrix_vector_product (u, r);
      if (alpha != 1.0)
#ifdef CSR_CALC_LIN_COMB_PARALLEL
  #pragma omp parallel for
#endif //CSR_CALC_LIN_COMB_PARALLEL
        for (i = 0; i < n; ++i)
          r[i] *= alpha;
    }
//  OMP_TIME_MEASURE_END (csr_calc_lin_comb_timer);
  return r_code;
}

/* ouput csr matrix to file
*  zero_symbol - symbol of zero element
*  format %2.2lf
*/
void
bubsort (double *vals, int *cols, int size)
{
 int i, j;
 double tmpv;
 int tmpc;
 for (i = 0; i < size - 1; i++)
   for (j = i + 1; j < size; j++)
     if (cols[i] > cols[j])
       {
       tmpc = cols[i];
       cols[i] = cols[j];
       cols[j] = tmpc;

       tmpv = vals[i];
       vals[i] = vals[j];
       vals[j] = tmpv;
       }
 }

int
csr_matrix::write_matrix_to_file (const char *file_name, const char *zero_symbol)
{
  FILE *fp = 0;
  int i, j, j1, j2;

  if (n_block_size > 1)
  {
    printf ("N_BLOCK_SIZE > 1!\n");
    return -1;
  }

  CH_PTR (file_name);
  CH_PTR (cols_ind);
  CH_PTR (rows_ptr);
  CH_PTR (values);

  fp = fopen (file_name, "w");
  if (!fp)
    //TODO: write error message
    return -1;

  fprintf (fp, "N_ROWS\tN_COLS\tN_NON_ZEROS\n");
  fprintf (fp, "%d\t%d\t%d\n\n", n_rows, n_cols, get_n_non_zeros ());
  fflush(fp);

  int *tmp=0;
  double* tmpv=0;
  int tmp_size = 0, r_code, flag;

  for (i = 0; i < n_cols; i++)
    fprintf (fp, "---------");
  fprintf (fp, "\n");

  fprintf (fp, "  \t%d\t", 0);
  for (i = 1; i < n_cols; i++)
    fprintf (fp, "%d\t", i);

  fprintf (fp, "\n");
  for (i = 0; i < n_cols; i++)
    fprintf (fp, "---------");
  fprintf (fp, "\n");

  for (i = 0; i < n_rows; i++)
    {
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];

      if (tmp_size < j2-j1)
        {
        tmp_size = j2-j1;
        FI_INT_ARRAY_REALLOCATOR (tmp, tmp_size, r_code);
        FI_DOUBLE_ARRAY_REALLOCATOR (tmpv, tmp_size, r_code);
        }

      for (j = j1; j < j2; j++)
         {
         tmp[j-j1] = cols_ind[j];
         tmpv[j-j1] = values[j];
         }

      bubsort (tmpv,tmp,j2-j1);

      fprintf (fp, "%d  |\t", i);

      for (int jj = 0; jj < n_cols; jj++)
      {
      flag = 0;
      for (j = j1; j < j2; j++)
        if(jj == tmp[j-j1])
          {
          fprintf (fp, "%.2lf\t", tmpv[j-j1]);
          flag = 1;
          }

      if(!flag) fprintf (fp, "%s\t", zero_symbol);

      }
      fprintf (fp, "\n");
    }

  fclose (fp);
  FI_FREE_ARRAY (tmp);
  FI_FREE_ARRAY (tmpv);
  return 0;
}

/*!
 * \brief write matrix to file with given name in ij format: ROW COLUMN VALUE
 *
 * \param file_name -- file to write
 *
 * \return 0 if success
 */
int
csr_matrix::write_matrix_to_file_ij (const char *file_name)
{
  FILE *fp = 0;
  int i, j, j1, j2, counter;

  if (n_block_size > 1)
  {
    printf ("N_BLOCK_SIZE > 1!\n");
    return -1;
  }

  CH_PTR (file_name);
  CH_PTR (cols_ind);
  CH_PTR (rows_ptr);
  CH_PTR (values);

  fp = fopen (file_name, "w");
  if (!fp)
    //TODO: write error message
    return -1;
 // fprintf (fp, "// N_ROWS\tN_COLS\tN_NON_ZEROS\tN_BLOCK_SIZE\n");
 // fprintf (fp, "%d\t%d\t%d\t%d\n", n_rows, n_cols, get_n_non_zeros (), n_block_size);

  for (i = 0, counter = 0; i < n_rows; ++i)
    {
      j1 = rows_ptr[i];
      j2 = rows_ptr[i + 1];
          for (j = j1; j < j2; ++j, ++counter)
            {
              fprintf (fp, "%d\t%d\t%.15lf\n", i + 1, cols_ind[j] + 1, values[j]);
            }
      fprintf (fp, "\n");
    }
  if (counter != get_n_non_zeros ())
    return -1;
  fclose (fp);
  return 0;
}


/*!
 * check for correctness structure of matrix(rows_ptr and cols_ind)
 */
int
csr_matrix::check_matrix ()
{
  int i, j, j1, j2, b_sqr, cl = 0, *columns = 0;
  int i_err = 0, jj = 0, j_last = 0, r_code = 0;

  CH_PTR (cols_ind);
  CH_PTR (rows_ptr);
  CH_PTR (values);

  b_sqr = n_block_size * n_block_size;
  int n_non_zeros = get_n_non_zeros ();
  printf ("check_matrix: n_rows=%d, n_cols=%d, nnz=%d, nb=%d\n", n_rows, n_cols, n_non_zeros, n_block_size);

  FI_INT_ARRAY_ALLOCATOR (columns, n_cols, r_code);
  memset (columns, -1, sizeof (int) * n_cols);

  if (rows_ptr[0] != 0)
    r_code = -1;

  for (i = 0; i < n_rows && !r_code; i++)
  {
    j1 = rows_ptr[i];
    j2 = rows_ptr[i + 1];

    if (j2 < j1)
    {
      i_err = i;
      r_code = -2;
      break;
    }

    if (j2 - j1 > n_cols)
    {
      i_err = i;
      r_code = -3;
      break;
    }

    if (cols_ind[j1] != i) //diagonal is first
    {
      i_err = i;
      r_code = -6;
      break;
    }

    j_last = jj;
    for (j = j1; j < j2; ++j)
    {
      cl = cols_ind[j];
      if (cl < 0 || cl >= n_cols)
      {
        i_err = i;
        r_code = -4;
        break;
      }
      if (columns[cl] >= j_last)
      {
        i_err = i;
        r_code = -5;
        break;
      }
      columns[cl] = jj++;
    }
  }

  if (!r_code)
    printf ("check_matrix: OK\n");
  else if (r_code == -1)
    printf ("rows_ptr[0] must be = 0!\n");
  else if (r_code == -2)
    printf ("in row %d: rows[i]>rows[i+1]! (%d,%d)\n", i_err, rows_ptr[i_err], rows_ptr[i_err+1]);
  else if (r_code == -3)
    printf ("in row %d: rows[i] - rows[i+1] > n_cols! (%d,%d)\n", i_err, rows_ptr[i_err], rows_ptr[i_err+1]);
  else if (r_code == -4)
  {
    printf ("in row %d: cols_ind must be in 0 <= .. < n_cols, but %d!\n", i_err, cl);
/*
    j1 = rows_ptr[i_err];
    j2 = rows_ptr[i_err + 1];
    for (j = j1; j < j2; ++j)
    printf ("%d\t", cols_ind[j]);
    printf ("\n");
*/
  }
  else if (r_code == -5)
    printf ("in row %d: duplicate cols (%d)!\n", i_err, cl);
  else if (r_code == -6)
    printf ("in row %d: diagonal is not first!\n", i_err);

  fflush(stdout);
  FI_FREE_ARRAY (columns);

  return r_code;
}
