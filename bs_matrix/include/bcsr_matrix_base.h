#ifndef __BCSR_MATRIX_BASE
#define __BCSR_MATRIX_BASE

#include "matrix_base.h"

namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   * <float, int>
   * <float, long>
   * <double, int>
   * <double, long>
   */
  template <class fp_vector_type, class i_vector_type>
  class BS_API_PLUGIN bcsr_matrix_base: public matrix_base<fp_vector_type, i_vector_type>
    {

    public:
      typedef typename fp_vector_type::value_type   fp_type_t;
      typedef typename i_vector_type::value_type    i_type_t;
      typedef fp_vector_type                        vec_fp_type;
      typedef fp_vector_type                        fp_vector_t;
      typedef i_vector_type                         i_vector_t;

      typedef fp_vector_type                        item_array_t;
      typedef i_vector_type                         index_array_t;
      typedef typename item_array_t::value_type     item_t;
      typedef typename index_array_t::value_type    index_t;

      typedef typename item_array_t::template array <float>::type   float_array_t;
      typedef typename item_array_t::template array <double>::type  double_array_t;

      typedef typename item_array_t::iterator       vec_iterator;
      typedef typename item_array_t::const_iterator vec_citerator;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL_T(bcsr_matrix_base);


    public:
      //! destructor
      virtual ~bcsr_matrix_base ()
      {
      }
      // initialize memory
      virtual int init (const this_t & /*matrix*/)
        {
          throw bs_exception ("init", "PURE CALL");
          return 0;
        }

      // initialize memory
      virtual int init (const i_type /*new_n_rows*/, const i_type /*new_n_cols*/, const i_type /*new_n_blok_size*/,
                         const i_type /*new_n_non_zeros*/)
        {
          throw bs_exception ("init", "PURE CALL");
          return 0;
        }

      // initialize memory
      virtual int init_struct (const i_type /*new_n_rows*/, const i_type /*new_n_cols*/, const i_type /*new_n_non_zeros*/)
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }


      // initialize indexes of main diagonal
      virtual int init_diag_ind (int diag_may_not_exist = 1, int row_offset = 0)
        {
          throw bs_exception ("init_diag_ind", "PURE CALL");
          return diag_may_not_exist + row_offset;
        }

      /**
       * \brief Copy matrix to current
       *
       * \param[in] matrix Source matrix
       * \return 0 if success
       */
      virtual int copy (this_t & /*matrix*/)
        {
          throw bs_exception ("copy", "PURE CALL");
          return 0;
        }

      // initialize only private arrays
      virtual int alloc_rows_ptr (const i_type /*new_n_rows*/, index_t diag_ind_default_value = 0)
        {
          throw bs_exception ("alloc_rows_ptr", "PURE CALL");
          return (int)diag_ind_default_value;
        }

      // allocate memory for column indexes and set it with default values
      virtual int alloc_cols_ind (const i_type /*new_n_non_zeros*/)
        {
          throw bs_exception ("alloc_cols_ind", "PURE CALL");
          return 0;
        }

      // allocate memory for values and set it with default values
      virtual int alloc_values (const i_type /*new_n_non_zeros*/)
        {
          throw bs_exception ("alloc_values", "PURE CALL");
          return 0;
        }

      //
      virtual int alloc_cols_ind_and_values (const i_type /*new_n_non_zeros*/)
        {
          throw bs_exception ("alloc_cols_ind_and_values", "PURE CALL");
          return 0;
        }

      //
      virtual int alloc_main_diagonal ()
        {
          throw bs_exception ("alloc_main_diagonal", "PURE CALL");
          return 0;
        }


      virtual int rand_init (const i_type /*new_n_rows*/, const i_type /*new_n_cols*/,
                     const fp_type /*rand_value_dispersion*/, const i_type /*elems_in_row*/)
        {
          throw bs_exception ("rand_init", "PURE CALL");
          return 0;
        }

      virtual int rand_init_symm (const i_type /*nx*/, const i_type /*ny*/, const i_type /*nz*/,
                          const fp_type /*rand_value_dispersion*/)
        {
          throw bs_exception ("rand_init_symm", "PURE CALL");
          return 0;
        }

      virtual int gen_2d_laplas (const i_type /*n*/)
        {
          throw bs_exception ("gen_2d_laplas", "PURE CALL");
          return 0;
        }
      // build B = A^T from given csr matrix, return 0 if success,
      // rows number and offset can be set manually with new_n_rows and rows_offset
      virtual int build_transpose (const this_t & /*matrix*/, const i_type rows_offset = 0,
                           const i_type cols_offset = 0, const i_type new_n_rows = 0)
        {
          throw bs_exception ("build_transpose", "PURE CALL");
          return (int)(rows_offset + cols_offset + new_n_rows);
        }

      // build B = A^T from given csr matrix using matrix structure only, without values, return 0 if success,
      // rows number and offset can be set manually with new_n_rows and rows_offset
      virtual int build_transpose_struct (const this_t &matrix,
                                  const i_type rows_offset = 0,
                                  const i_type cols_offset = 0,
                                  const i_type new_n_rows = 0);
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

      //! return number of nonzeros elements
      virtual i_type get_n_non_zeros () const
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

      //! check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check ();
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

#ifdef _HDF5_OUT
      // calculate matrix vector product using hdf5 dump, v -- input vector, r -- output vector,
      // hyperslab_size -- size of hyperslab to be read from file
      // r += A * v
      // return 0 if success
      virtual int matrix_vector_product_hdf5 (const fp_type *v, fp_type *r, i_type hyperslab_size);
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

      // calculate matrix vector product using hdf5 dump (multithreaded), v -- input vector, r -- output vector,
      // hyperslab_size -- size of hyperslab to be read from file
      // r += A * v
      // return 0 if success
      virtual int matrix_vector_product_hdf5_mt (const fp_type *v, fp_type *r, i_type hyperslab_size);
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }
#endif // #ifdef _HDF5_OUT

      ////////////////////////
      //  IO methods
      ////////////////////////

      // write matrix to file with given name with or without sorting cols
      virtual int ascii_write_in_csr_format (const std::string &file_name) const;
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

      // print matrix to ascii file in IJ format
      virtual int ascii_write_in_ij_format (const std::string &file_name) const;
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

      // read matrix from file with given name (csr format)
      virtual int ascii_read_from_csr_format (const std::string &file_name);
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

#ifdef _HDF5_OUT
      // dump cols_ind and values to hdf5 file
      virtual int dump_to_hdf5 ();
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

      // read hyperslabs of cols_ind and values from hdf5 dump
      virtual int read_from_hdf5(i_type hyperslab_size,
                         i_vector_type cols_ind_hyperslab[],
                         fp_vector_type values_hyperslab[],
                         boost::mutex mutex[],
                         boost::condition cond[],
                         bool data_ready[]);
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }
#endif //_HDF5 // #ifdef _HDF5

      //! merge two csr matrices into current
      virtual int merge(const this_t *lhs, const this_t *rhs, int key, index_array_t &cfl_vector);
        {
          throw bs_exception ("init_struct", "PURE CALL");
          return 0;
        }

      virtual index_array_t &
      get_rows_ptr ()
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual const index_array_t &
      get_rows_ptr () const
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual index_array_t &
      get_cols_ind ()
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual const index_array_t &
      get_cols_ind () const
      {
        static index_array_t dummy;
        return dummy;
      }

      virtual item_array_t &
      get_values ()
      {
        static item_array_t dummy;
        return dummy;
      }

      virtual const item_array_t &
      get_values () const
      {
        static item_array_t dummy;
        return dummy;
      }
    public:
      ////////////////////////
      // DEBUG METHODS
      ///////////////////////
      int print_statistic ();
    };

}//namespace blue_sky
#endif //__BCSR_MATRIX_BASE

