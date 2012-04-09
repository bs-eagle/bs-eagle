#ifndef __MAT_BASE_H
#define __MAT_BASE_H

/*! 
 * \file mat_base.h
 * \brief base class for matrix representation
 * \author Borschuk Oleg
 * \date 2006-07-26
 */

enum
{
  MATRIX_TYPE_BANDED = 1,
  MATRIX_TYPE_CSR,
  MATRIX_TYPE_2_IN_1,
  MATRIX_TYPE_MPI_CSR_SIMPLE,
  MATRIX_TYPE_MPI_CSR_DIAGOFFD,
  MATRIX_TYPE_MPI_2_IN_1
};

class mat_base
{
  //-----------------------------------------
  //  METHODS
  //-----------------------------------------
  public:
    //! constructor
    mat_base () 
      {
        type = -1;
        n_block_size = 1;
        n_rows = 0;
        n_cols = 0;
        is_square = 0;
        n_row_size = 1;
      }
    //! destructor 
    virtual ~mat_base () 
      {
        type = -1;
        n_block_size = 1;
        n_rows = 0;
        n_cols = 0;
        is_square = 0;
        n_row_size = 1;
      }

    // calculate matrix vector product, v -- input vector, r -- output vector
    // r += A * v
    // return 0 if success
    virtual int matrix_vector_product (double *v, double *r) = 0;

  private:
  //-----------------------------------------
  //  VARIABLES
  //-----------------------------------------
  public:
    int type;
    int n_block_size;
    int n_rows;
    int n_cols;
    int is_square;
    int n_row_size;     //! number of rows in each block
    
#ifdef _MPI
    int comm;          //!< global MPI communicator
    
    int *row_starts;        //!< array [n_procs + 1] of row start indexes on each proc + total row count
    int row_start;          //!< row offset of current proc row_starts[proc_num]
       
    int local_n_rows;       //!< local row count
    int global_n_rows;      //!< global row count
    
    // used for MATRIX_TYPE_MPI_CSR_DIAGOFFD type
    int *col_starts;        //!< array [n_procs + 1] of col start indexes on each proc + total col count
#endif //_MPI

  private:
};


#endif //__MAT_BASE_H
