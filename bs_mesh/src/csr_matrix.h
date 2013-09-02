#ifndef __CSR_MATRIX__H
#define __CSR_MATRIX__H
/*!
* \file csr_matrix.h
* \brief declaration of compressed sparce row matrix without diagonal
* \author Borschuk Oleg
* \date 2006-06-16
*/

#include "mat_base.h"


//-----------------------------------------------------------
//-----------------------------------------------------------
class csr_matrix : public mat_base
{
  //-----------------------------------------
  //  METHODS
  //-----------------------------------------
  public:
    // default constructor
    csr_matrix ();

    // default destructor
    ~csr_matrix ();

    // initialize memory
    int init (csr_matrix *matrix);
    void setup_all (const int new_n_rows, const int new_n_block_size, double *new_values, int *new_rows_ptr, int *new_cols_ind);

    // initialize memory
    int init (const int new_n_rows, const int new_n_cols, const int new_n_blok_size,
              const int new_n_non_zeros);

    // initialize memory
    int init_struct (const int new_n_rows, const int new_n_cols, const int new_n_non_zeros);


    // initialize only private arrays
    int alloc_rows_ptr (const int new_n_rows);
    int alloc_cols_ind (const int new_n_non_zeros);
    int alloc_values (const int new_n_non_zeros);
    int alloc_cols_ind_and_values (const int new_n_non_zeros);
    int alloc_diag_ind ();

    // initialize diagonal indexes
    int init_diag_ind (int diag_may_not_exist = 1, int row_offset = 0);
    
    // initialize next to diagonal indexes (correct only with sorted matrix)
    int init_next_to_diag_ind (int row_offset = 0, int *col_map = 0);

    int rand_init (const int new_n_rows, const int new_n_cols,
                   const double rand_value_dispersion, int elems_in_row);
    int rand_init_symm (const int nx, const int ny, const int nz, const int nb,
                        const double rand_value_dispersion);
    int gen_2d_laplas (const int n);
    // calculate matrix vector product, v -- input vector, r -- output vector
    // r += A * v
    int matrix_vector_product (double *v, double *r);

    // transpose matrix vector product
    // r += A^T * v
    int matrix_vector_product_t (double *v, double *r);

    // calculate linear combination r = alpha * Au + beta * v
    int calc_lin_comb (const double alpha, const double beta,
                       double *u, double *v, double *r);

    // build B = A^T from given csr matrix, return 0 if success,
    // rows number and offset can be set manually with new_n_rows and rows_offset
    int build_transpose (csr_matrix *matrix, int rows_offset = 0, int cols_offset = 0, int new_n_rows = 0);

    // build B = A^T from given csr matrix using matrix structure only, without values, return 0 if success,
    // rows number and offset can be set manually with new_n_rows and rows_offset
    int build_transpose_struct (csr_matrix *matrix, int rows_offset = 0, int cols_offset = 0, int new_n_rows = 0);

    //! return values array
    double *get_values () const {return values;}
    //! return rows_ptr array
    int *get_rows_ptr () const {return rows_ptr;}
    //! return cols_ind array
    int *get_cols_ind () const {return cols_ind;}
    //! return diag_ind array
    int *get_diag_ind () const {return diag_ind;}
    //void set_data (double *new_values, int *new_rows_ptr, int *new_cols_ind);

    //! return number of nonzeros elements
    int get_n_non_zeros ();

    //! check for correctness the structure of matrix (rows_ptr and cols_ind)
    int check_matrix ();

    //! read matrix in ascii format
    int read_matrix_in_ascii_format (char *file);


    ////////////////////////
    //  IO methods
    ////////////////////////
    // write matrix to file with given name with or without sorting cols
    int write_matrix_to_file (const char *file_name, int sort_cols = 0);
    int write_matrix_to_file (const char *file_name, const char *zero_symbol);
    int mpi_write_matrix_to_file (const char *file_name, int sort_cols = 0);
    int write_matrix_to_file_ij (const char *file_name);

    // read matrix from file with given name
    int read_matrix_from_file (const char *file_name);

    ////////////////////////
    // DEBUG METHODS
    ///////////////////////
    int print_statistic ();

  //-----------------------------------------
  //  VARIABLES
  //-----------------------------------------
  public:

  private:
    int n_memory_n_rows;                //!< amount of memory allocated for rows_ptr and diag_ind
    int n_memory_cols_ind;              //!< amount of memory allocated for cols_ind
    int n_memory_values;                //!< amount of memory allocated for values

    int own_memory;                     //!< specifies delete or not data storages in destructor

    // matrix storage
    double *values;                     //!< storage for block matrix values (compressed row base 0 format)

    // matrix indexes
    int *rows_ptr;                      //!< start row indexes for block matrix
    int *cols_ind;                      //!< column indexes for all blocks in matrix
    int *diag_ind;                      //!< indexes of diagonals
};

void
bubsort (double *vals, int *cols, int size);
#endif //__CSR_MATRIX__H
