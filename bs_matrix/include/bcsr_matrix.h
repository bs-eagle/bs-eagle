/**
 * @file bcsr_matrix.h
 * @brief
 * @author Borschuk Oleg
 * @date 2008-03-25
 */

#ifndef BCSR_MATRIX_H_
#define BCSR_MATRIX_H_

#include "bs_object_base.h"

#include "matrix_base.h"
#include "matrix_macroses.h"
#include "setup_preconditioner.h"
#include "shared_vector.h"

#include "vector_assign.h"

namespace blue_sky
{



  template <class fp_vector_t, class i_vector_t>
  class BS_API_PLUGIN bcsr_matrix : public matrix_base<fp_vector_t, i_vector_t>//, public setup_preconditioner<matrix_base<fp_vector_t, i_vector_t> >
    {
      template <class,class >
      friend class py_bcsr_matrix;

    public:

      typedef matrix_base <fp_vector_t, i_vector_t> base_t;
      typedef base_t                                matrix_t;

      using matrix_t::n_rows;
      using matrix_t::n_cols;
      using matrix_t::n_block_size;

      //-----------------------------------------
      //  TYPE DECL
      //-----------------------------------------
    public:
      typedef fp_vector_t                           fp_vector_type;
      typedef i_vector_t                            i_vector_type;
      typedef typename fp_vector_t::value_type      fp_type;
      typedef typename i_vector_t::value_type       i_type;

      typedef fp_vector_t                           item_array_t;

      typedef i_vector_t                            index_array_t;
      typedef typename item_array_t::value_type     item_t;
      typedef typename index_array_t::value_type    index_t;

      typedef bcsr_matrix<fp_vector_t, i_vector_t>  this_t;
      typedef matrix_t                              matrix_base_t;

      typedef shared_vector <float>                 float_array_t;
      typedef shared_vector <double>                double_array_t;

      typedef smart_ptr <matrix_t, true>            sp_matrix_t;
      typedef smart_ptr<this_t, true>               sp_bcsr_matrix_t;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL(bcsr_matrix);

    // TODO: HACK
    public:
      struct hack_bcsr_matrix {};
      bcsr_matrix (const hack_bcsr_matrix &);

    public:
      // default constructor
      //bcsr_matrix ();

      // default destructor
      virtual ~bcsr_matrix ();

      // initialize memory
      int init (const this_t &matrix);

      // initialize memory
      int init (const i_type new_n_rows, const i_type new_n_cols, const i_type new_n_blok_size,
                const i_type new_n_non_zeros);

      // initialize memory
      int init_struct (const i_type new_n_rows, const i_type new_n_cols, const i_type new_n_non_zeros);
      
      // initialize memory
      int init_values (const i_type new_n_block_size);

      // initialize indexes of main diagonal
      int init_diag_ind (int diag_may_not_exist = 1, int row_offset = 0);

      /**
       * \brief Copy matrix to current
       *
       * \param[in] matrix Source matrix
       * \return 0 if success
       */
      int copy (this_t &matrix);

      // initialize only private arrays
      int alloc_rows_ptr (const i_type new_n_rows, index_t diag_ind_default_value = 0);

      // allocate memory for column indexes and set it with default values
      int alloc_cols_ind (const i_type new_n_non_zeros);

      // allocate memory for values and set it with default values
      int alloc_values (const i_type new_n_non_zeros);

      //
      int alloc_cols_ind_and_values (const i_type new_n_non_zeros);

      //
      int alloc_main_diagonal ();


      int rand_init (const i_type new_n_rows, const i_type new_n_cols,
                     const fp_type rand_value_dispersion, const i_type elems_in_row);

      int rand_init_symm (const i_type nx, const i_type ny, const i_type nz,
                          const fp_type rand_value_dispersion);

      int gen_2d_laplas (const i_type n);
      // build B = A^T from given csr matrix, return 0 if success,
      // rows number and offset can be set manually with new_n_rows and rows_offset
      int build_transpose (const this_t &matrix, const i_type rows_offset = 0,
                           const i_type cols_offset = 0, const i_type new_n_rows = 0);

      // build B = A^T from given csr matrix using matrix structure only, without values, return 0 if success,
      // rows number and offset can be set manually with new_n_rows and rows_offset
      int build_transpose_struct (const this_t &matrix,
                                  const i_type rows_offset = 0,
                                  const i_type cols_offset = 0,
                                  const i_type new_n_rows = 0);

      //! return values array
      fp_vector_type &get_values ()
      {
        return values;
      }

      const fp_vector_type &get_values () const
        {
          return values;
        }

      //! return rows_ptr array
      i_vector_type &get_rows_ptr ()
      {
        return rows_ptr;
      }

      //! return rows_ptr array
      const i_vector_type &get_rows_ptr () const
      {
        return rows_ptr;
      }

      //! return cols_ind array
      i_vector_type &get_cols_ind ()
      {
        return cols_ind;
      }

      //! return cols_ind array
      const i_vector_type &get_cols_ind () const
        {
          return cols_ind;
        }

      //! return diag_ind array
      i_vector_type &get_diag_ind ()
      {
        return diag_ind;
      }

      //! return diag_ind array for read only
      const i_vector_type &get_diag_ind () const
        {
          return diag_ind;
        }


      //! return number of nonzeros elements
      i_type get_n_non_zeros () const
        {
          BS_ASSERT (base_t::n_rows < (int)rows_ptr.size ()) (base_t::n_rows) ((int)rows_ptr.size ());
          return rows_ptr.empty () ? 0 : rows_ptr.back ();
        }

      //! check for correctness the structure of matrix (rows_ptr and cols_ind)
      int internal_check ();

#ifdef _HDF5_OUT
      // calculate matrix vector product using hdf5 dump, v -- input vector, r -- output vector,
      // hyperslab_size -- size of hyperslab to be read from file
      // r += A * v
      // return 0 if success
      virtual int matrix_vector_product_hdf5 (const fp_type *v, fp_type *r, i_type hyperslab_size);

      // calculate matrix vector product using hdf5 dump (multithreaded), v -- input vector, r -- output vector,
      // hyperslab_size -- size of hyperslab to be read from file
      // r += A * v
      // return 0 if success
      virtual int matrix_vector_product_hdf5_mt (const fp_type *v, fp_type *r, i_type hyperslab_size);
#endif // #ifdef _HDF5_OUT

      // calculate matrix vector product, v -- input vector, r -- output vector
      // r += A * v
      // return 0 if success
      virtual int matrix_vector_product (const double_array_t &v, float_array_t &r) const;
      virtual int matrix_vector_product (const float_array_t &v, float_array_t &r) const;
      virtual int matrix_vector_product (const double_array_t &v, double_array_t &r) const;

      // calculate matrix^t vector product, v -- input vector, r -- output vector
      // r += A^t * v
      // return 0 if success
      virtual int matrix_vector_product_t (const double_array_t &v, float_array_t &r);
      virtual int matrix_vector_product_t (const float_array_t &v, float_array_t &r);
      virtual int matrix_vector_product_t (const double_array_t &v, double_array_t &r);

      // calculate linear combination r = alpha * Au + beta * v
      virtual int calc_lin_comb (item_t                 alpha,
                                 item_t                 beta,
                                 const double_array_t   &u,
                                 const float_array_t    &v,
                                 double_array_t         &r) const;
      virtual int calc_lin_comb (item_t                 alpha,
                                 item_t                 beta,
                                 const float_array_t    &u,
                                 const double_array_t   &v,
                                 float_array_t          &r) const;
      virtual int calc_lin_comb (item_t                 alpha,
                                 item_t                 beta,
                                 const float_array_t    &u,
                                 const float_array_t    &v,
                                 float_array_t          &r) const;
      virtual int calc_lin_comb (item_t                 alpha,
                                 item_t                 beta,
                                 const double_array_t   &u,
                                 const double_array_t   &v,
                                 double_array_t         &r) const;

      // return total amount of allocated memory
      virtual i_type get_allocated_memory_in_bytes ();

      ////////////////////////
      //  IO methods
      ////////////////////////

      // write matrix to file with given name with or without sorting cols
      int ascii_write_in_csr_format (const std::string &file_name) const;

      // print matrix to ascii file in IJ format
      int ascii_write_in_ij_format (const std::string &file_name) const;

      // read matrix from file with given name (csr format)
      int ascii_read_from_csr_format (const std::string &file_name);

#ifdef _HDF5_OUT
      // dump cols_ind and values to hdf5 file
      int dump_to_hdf5 ();

      // read hyperslabs of cols_ind and values from hdf5 dump
      int read_from_hdf5(i_type hyperslab_size,
                         i_vector_type cols_ind_hyperslab[],
                         fp_vector_type values_hyperslab[],
                         boost::mutex mutex[],
                         boost::condition cond[],
                         bool data_ready[]);
#endif //_HDF5 // #ifdef _HDF5

      //! merge two csr matrices into current
      int merge(const this_t *lhs, const this_t *rhs, int key, const item_array_t &cfl_vector);

      virtual void
      init_vector (float_array_t &v)
      {
        v.assign (this->n_rows * this->n_block_size, 0);
      }
      virtual void
      init_vector (double_array_t &v)
      {
        v.assign (this->n_rows * this->n_block_size, 0);
      }
      virtual void
      init_vector (i_vector_t &v)
      {
        v.assign (this->n_rows * this->n_block_size, 0);
      }

    public:
      ////////////////////////
      // DEBUG METHODS
      ///////////////////////
      int print_statistic ();

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:

    private:
      fp_vector_type values;              //!< matrix values (nnz * n_block_size * n_block_size)

      i_vector_type diag_ind;             //!< diagonal indexies
      i_vector_type rows_ptr;             //!< row start indexies (n_rows + 1), rows_ptr[0] = 0
      i_vector_type cols_ind;             //!< column indexies for all blocks (nnz)

    };



  namespace python
  {
    void py_export_base_matrices();
  }

} //ns bs

#endif//BCSR_MATRIX_H_

