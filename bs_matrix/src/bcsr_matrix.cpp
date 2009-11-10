/**
 * @file bcsr_matrix.cpp
 * @brief
 * @author Borschuk Oleg
 * @date 2008-03-25
 */
#include "bs_matrix_stdafx.h"

#include "bcsr_matrix.h"
#include "py_bcsr_matrix.h"
#include "locale_keeper.h"
#include "bos_report.h"

#include "merge_matrices.h"
#include "bcsr_matrix_vector_product.h"
#include "bcsr_calc_lin_comb.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
  {
#ifdef _HDF5_OUT
  const char hdf5_dump_path[] = "d:\\temp\\temp.hdf5";
#endif // #ifdef _HDF5_OUT

  //! default constructor
  template<class fp_vector_t, class i_vector_t>
  bcsr_matrix<fp_vector_t, i_vector_t>::bcsr_matrix(bs_type_ctor_param)
  {
  }

  // TODO: HACK
  template <class fp_vector_t, class i_vector_t>
  bcsr_matrix <fp_vector_t, i_vector_t>::bcsr_matrix (const hack_bcsr_matrix &)
  {
  }

  template< class fp_vector_t,class i_vector_t >
  bcsr_matrix<fp_vector_t, i_vector_t>::bcsr_matrix(const this_t &src)
        : bs_refcounter (), matrix_base<fp_vector_t, i_vector_t> ()
  {
    *this = src;
  }

  //! default destructor
  template<class fp_vector_t, class i_vector_t>
  bcsr_matrix<fp_vector_t, i_vector_t>::~bcsr_matrix ()
  {
  }

  /*!
    \brief initialize matrix
    \param new_n_rows -- number of matrix rows (in blocks)
    \param new_n_cols -- number of matrix columns (in blocks)
    \param new_n_blok_size -- size of block
    \param new_n_non_zeros -- number of non zeros element in matrix (in blocks)
    \return 0 if success
  */
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::init (const i_type new_n_rows,
      const i_type new_n_cols,
      const i_type new_n_block_size,
      const i_type new_n_non_zeros)
  {
    if (new_n_block_size < 1 || new_n_cols < 1)
      return -1;

    this->n_block_size = new_n_block_size;
    this->n_cols = new_n_cols;

    if (alloc_rows_ptr (new_n_rows)
        || alloc_cols_ind_and_values (new_n_non_zeros))
      {
        return -2;
      }

    return 0;
  }

  /**
   * \brief Copy matrix to current
   *
   * \param[in] matrix Source matrix
   * \return 0 if success
   */
  template <class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::copy (this_t &matrix)
  {
    BS_ASSERT (base_t::n_block_size);
    BS_ASSERT (base_t::n_cols);
    BS_ASSERT (base_t::n_rows);

    if (!base_t::n_block_size)
      return -1;
    if (!base_t::n_cols)
      return -2;
    if (!base_t::n_rows)
      return -3;

    rows_ptr.assign (matrix.get_rows_ptr ().begin (), matrix.get_rows_ptr ().end ());
    cols_ind.assign (matrix.get_cols_ind ().begin (), matrix.get_cols_ind ().end ());
    values.assign   (matrix.get_values ().begin (), matrix.get_values ().end ());
    assign (diag_ind, 0);

    BS_ASSERT (rows_ptr.size () == (size_t)base_t::n_rows);
    BS_ASSERT (cols_ind.size () == (size_t)base_t::n_cols);
    BS_ASSERT (values.size ()   == (size_t)base_t::n_cols * base_t::n_block_size);

    return 0;
  }

  /*!
   */
  template <class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::alloc_rows_ptr (const i_type new_n_rows, index_t diag_ind_value)
  {
    if (new_n_rows < 1)
      return -1;

    n_rows = new_n_rows;

    //rows_ptr.assign (n_rows + 1, 0);			// miryanov at 11.09.2008 - jacobian_matrix->init
    // new size is n_rows + 1
    rows_ptr.resize (n_rows + 1, rows_ptr.size () ? rows_ptr.back () : 0);
    assign (diag_ind, n_rows, diag_ind_value);

    return 0;
  }

  /*!
   */
  template <class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::alloc_cols_ind (const i_type new_n_non_zeros)
  {
    int nnz = new_n_non_zeros;

    if (nnz < 1)
      nnz = 1;

    assign (cols_ind, nnz, -1);
    return 0;
  }

  /*!
   */
  template <class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::alloc_values (const i_type new_n_non_zeros)
  {
    int b_sqr = n_block_size * n_block_size;
    int nnz = new_n_non_zeros;

    if (nnz < 1)
      nnz = 1;

    assign (values, b_sqr * nnz, 0);

    return 0;
  }

  /*!
   */
  template <class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::alloc_cols_ind_and_values (const i_type new_n_non_zeros)
  {
    if (alloc_cols_ind (new_n_non_zeros) || alloc_values (new_n_non_zeros))
      return -2;

    return 0;
  }

  /*!
   * \brief initialize diagonal indexes
   *
   * \return 0 if success
   */
  template<class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix<fp_vector_t, i_vector_t>::init_diag_ind (int diag_may_not_exist, int row_offset)
  {
    if (rows_ptr.empty ())
      return 0;

    assign (diag_ind, n_rows, -1);

    for (index_t i = 0, cnt = (index_t)rows_ptr.size () - 1; i < cnt; ++i)
      {
        index_t j1 = rows_ptr[i];
        index_t j2 = rows_ptr[i + 1];

        for (index_t j = j1; j < j2; ++j)
          {
            if (cols_ind[j] == i + row_offset)
              {
                diag_ind[i] = j;
                break;
              }
          }

        if (diag_ind[i] == -1 && !diag_may_not_exist) // diagonal not found
          return -4;
      }

    return 0;
  }
  //! init by parent matrix
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::init (const this_t &matrix)
  {
    if (init (matrix.n_rows, matrix.n_cols, matrix.n_block_size, matrix.get_n_non_zeros ()))
      return -2;

    rows_ptr.assign (matrix.get_rows_ptr ().begin (), matrix.get_rows_ptr ().end ());
    cols_ind.assign (matrix.get_cols_ind ().begin (), matrix.get_cols_ind ().end ());
    values.assign   (matrix.get_values ().begin (), matrix.get_values ().end ());
    return 0;
  }

  /*!
    \brief initialize matrix struct
    \param new_n_rows -- number of matrix rows (in blocks)
    \param new_n_cols -- number of matrix columns (in blocks)
    \param new_n_non_zeros -- number of non zeros element in matrix (in blocks)
    \return 0 if success
  */
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::init_struct (const i_type new_n_rows, const i_type new_n_cols, const i_type new_n_non_zeros)
  {
    if (new_n_cols < 1)
      return -1;

    n_cols = new_n_cols;
    if (alloc_rows_ptr (new_n_rows) || alloc_cols_ind (new_n_non_zeros))
      return -2;

    return 0;
  }
  
  /*!
    \brief initialize matrix values
    \param new_n_block_size -- dimension of single block
    \return 0 if success
  */
  template<class fp_vector_t, class i_vector_t> int
    bcsr_matrix<fp_vector_t, i_vector_t>::init_values (const i_type new_n_block_size)
  {
    if (new_n_block_size < 1)
      return -1;

    n_block_size = new_n_block_size;
    if (alloc_values ((index_t)cols_ind.size()))
      return -2;

    return 0;
  }

#ifdef _HDF5_OUT
  /*!
   * \brief matrix vector product using hdf5 dump
   * \param v -- vector
   * \param r -- result vector
    * \param hyperslab_size
   * \return 0 if success
   */
  template<class fp_type, class i_type> int
  bcsr_matrix<fp_type, i_type>::matrix_vector_product_hdf5 (const fp_type *v, fp_type *r, i_type hyperslab_size)
  {
    H5::H5File *file = new H5::H5File(hdf5_dump_path, H5F_ACC_RDONLY);
    H5::DataSet dataset_cols_ind = file->openDataSet("/cols_ind");
    H5::DataSet dataset_values = file->openDataSet("/values");

    i_type i, j, k;
    i_type j1, j2;
    i_type cl;
    i_type nnz = get_n_non_zeros();
    i_vector_type cols_ind_hyperslab;
    fp_vector_type values_hyperslab;
    fp_type *r_block, *m_block;
    const fp_type *v_block;
    int b_sqr = n_block_size * n_block_size;

    OMP_TIME_MEASURE_START (csr_matrix_vector_product_timer);

    CH_PTR (v);
    CH_PTR (r);

    ////////////////////////////////////////////
    // loop through rows
    ////////////////////////////////////////////

#ifdef CSR_MATRIX_VECTOR_PRODUCT_PARALLEL
#pragma omp parallel for private (r_block, j1, j2, j, k, cl, m_block, v_block)
#endif //CSR_MATRIX_VECTOR_PRODUCT_PARALLEL

    for (i = 0; i < n_rows; ++i)
      {
        r_block = r + i * n_block_size;
        j1 = rows_ptr[i];
        j2 = rows_ptr[i + 1];

        for (j = j1; j < j2; j += hyperslab_size)
          {
            i_type real_hyperslab_size = min(hyperslab_size, j2 - j);

            // reading data from dump
            hsize_t dims_memory[] = {real_hyperslab_size};
            H5::DataSpace dataspaceMemory(1, dims_memory);
            hsize_t dims_file[] = {nnz};
            H5::DataSpace dataspaceFile(1, dims_file);
            hsize_t offset[] = {j};
            hsize_t count[] = {real_hyperslab_size};
            dataspaceFile.selectHyperslab(H5S_SELECT_SET, count, offset);
            cols_ind_hyperslab.resize(real_hyperslab_size);
            dataset_cols_ind.read(&cols_ind_hyperslab[0], H5::PredType::NATIVE_INT, dataspaceMemory, dataspaceFile);
            values_hyperslab.resize(real_hyperslab_size);
            dataset_values.read(&values_hyperslab[0], H5::PredType::NATIVE_FLOAT, dataspaceMemory, dataspaceFile);

            // calc
            for (k = 0; k < real_hyperslab_size; ++k)
              {
                cl = cols_ind_hyperslab[k];
                m_block = &values_hyperslab[k * b_sqr];
                v_block = v + cl * n_block_size;
                MV_PROD (n_block_size, m_block, v_block, r_block);
              }
          }
      }

    OMP_TIME_MEASURE_END (csr_matrix_vector_product_timer);
    return 0;
  }

  /*!
  * \brief matrix vector product using hdf5 dump (multithreaded)
  * \param v -- vector
  * \param r -- result vector
   * \param hyperslab_size
  * \return 0 if success
  */
  template<class fp_type, class i_type> int
  bcsr_matrix<fp_type, i_type>::matrix_vector_product_hdf5_mt (const fp_type *v, fp_type *r, i_type hyperslab_size)
  {
    i_type i, j, k, l;
    i_type j1, j2;
    i_type cl;
    fp_type *r_block, *m_block;
    const fp_type *v_block;
    int b_sqr = n_block_size * n_block_size;

    OMP_TIME_MEASURE_START (csr_matrix_vector_product_timer);

    CH_PTR (v);
    CH_PTR (r);

    i_vector_type cols_ind_hyperslab[2];
    fp_vector_type values_hyperslab[2];
    boost::mutex mutex[2];
    boost::condition cond[2];
    bool data_ready[2];
    data_ready[0] = false;
    data_ready[1] = false;

    // launch thread for reading from hdf5 dump
    boost::thread threadRead(boost::bind(&this_t::read_from_hdf5,
                                         boost::ref(*this),
                                         hyperslab_size,
                                         cols_ind_hyperslab,
                                         values_hyperslab,
                                         mutex,
                                         cond,
                                         data_ready));

    ////////////////////////////////////////////
    // loop through rows
    ////////////////////////////////////////////

#ifdef CSR_MATRIX_VECTOR_PRODUCT_PARALLEL
#pragma omp parallel for private (r_block, j1, j2, j, k, l, cl, m_block, v_block)
#endif //CSR_MATRIX_VECTOR_PRODUCT_PARALLEL
    for (i = 0, k = 0; i < n_rows; ++i)
      {
        r_block = r + i * n_block_size;
        j1 = rows_ptr[i];
        j2 = rows_ptr[i + 1];

        for (j = j1; j < j2; j += hyperslab_size)
          {
            {
              // waiting while current data region is blocked
              boost::mutex::scoped_lock lock(mutex[k]);
              while (data_ready[k] == false)
                cond[k].wait(lock);
              // calc using data from current data region
              for (l = 0; l < cols_ind_hyperslab[k].size(); ++l)
                {
                  cl = cols_ind_hyperslab[k][l];
                  m_block = &values_hyperslab[k][l * b_sqr];
                  v_block = v + cl * n_block_size;
                  MV_PROD (n_block_size, m_block, v_block, r_block);
                }
              data_ready[k] = false;
            }
            cond[k].notify_one();
            // switching to the next data region
            k = (k + 1) % 2;
          }
      }

    OMP_TIME_MEASURE_END (csr_matrix_vector_product_timer);
    return 0;
  }
#endif // #ifdef _HDF5_OUT

  /*!
    \brief matrix vector product
    \param v -- vector
    \param r -- result vector
    \return 0 if success
  */
  template<class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix<fp_vector_t, i_vector_t>::matrix_vector_product (const double_array_t &v, float_array_t &r) const
  {
    return bcsr_matrix_vector_product (this, v, r);
  }
  template<class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix<fp_vector_t, i_vector_t>::matrix_vector_product (const float_array_t &v, float_array_t &r) const
  {
    return bcsr_matrix_vector_product (this, v, r);
  }
  template<class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix<fp_vector_t, i_vector_t>::matrix_vector_product (const double_array_t &v, double_array_t &r) const
  {
    return bcsr_matrix_vector_product (this, v, r);
  }

  /*!
   * \brief calculate r += A^T * v
   *
   * \param v -- vector
   * \param r -- result vector
   *
   * \return 0 -- if success
   */
  template <class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix<fp_vector_t, i_vector_t>::matrix_vector_product_t (const double_array_t &v, float_array_t &r)
  {
    return bcsr_matrix_vector_product_t (this, v, r);
  }
  template <class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix<fp_vector_t, i_vector_t>::matrix_vector_product_t (const float_array_t &v, float_array_t &r)
  {
    return bcsr_matrix_vector_product_t (this, v, r);
  }
  template <class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix<fp_vector_t, i_vector_t>::matrix_vector_product_t (const double_array_t &v, double_array_t &r)
  {
    return bcsr_matrix_vector_product_t (this, v, r);
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
  template <class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix <fp_vector_t, i_vector_t>::calc_lin_comb (item_t alpha, item_t beta, const double_array_t &u, const float_array_t &v, double_array_t &r) const
  {
    return bcsr_calc_lin_comb (this, alpha, beta, u, v, r);
  }
  template <class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix <fp_vector_t, i_vector_t>::calc_lin_comb (item_t alpha, item_t beta, const float_array_t &u, const double_array_t &v, float_array_t &r) const
  {
    return bcsr_calc_lin_comb (this, alpha, beta, u, v, r);
  }
  template <class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix <fp_vector_t, i_vector_t>::calc_lin_comb (item_t alpha, item_t beta, const float_array_t &u, const float_array_t &v, float_array_t &r) const
  {
    return bcsr_calc_lin_comb (this, alpha, beta, u, v, r);
  }
  template <class fp_vector_t, class i_vector_t>
  int
  bcsr_matrix <fp_vector_t, i_vector_t>::calc_lin_comb (item_t alpha, item_t beta, const double_array_t &u, const double_array_t &v, double_array_t &r) const
  {
    return bcsr_calc_lin_comb (this, alpha, beta, u, v, r);
  }

  // return total amount of allocated memory
  template<class fp_vector_t, class i_vector_t> typename bcsr_matrix<fp_vector_t, i_vector_t>::i_type
  bcsr_matrix<fp_vector_t, i_vector_t>::get_allocated_memory_in_bytes ()
  {
    i_type n = 0;

    n += (i_type)sizeof (*this);
    n += (i_type)rows_ptr.capacity () * (i_type)sizeof (i_type);
    n += (i_type)cols_ind.capacity () * (i_type)sizeof (i_type);
    n += (i_type)diag_ind.capacity () * (i_type)sizeof (i_type);
    n += (i_type)values.capacity () * (i_type)sizeof (fp_type);

    return n;
  }


  /**
   * @brief build transpose matrix
   *
   * @param matrix        -- <INPUT> given base matrix for bilding transpose matrix
   * @param rows_offset   --
   * @param cols_offset   --
   * @param new_n_rows    --
   *
   * @return
   */
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::build_transpose (const this_t &matrix,
      const i_type rows_offset,
      const i_type cols_offset,
      const i_type new_n_rows)
  {
    i_type i,j;
    i_type row_ind, j1, j2;
    i_type nnz;
    i_type nr;

    nnz = matrix.get_n_non_zeros ();

    if (new_n_rows == 0)
      nr = matrix.n_cols;
    else
      nr = new_n_rows;


    init (nr, matrix.n_rows, matrix.n_block_size, nnz);

    // TODO: BUG:
    assign (rows_ptr, 0);
    // Count the number of entries in each column of matrix (row of transposed matrix) and fill the rows_ptr array.

    for (i = 0; i < nnz; i++)
      {
        ++rows_ptr[matrix.cols_ind[i] + 1 - rows_offset];
      }

    for (i = 2; i <= nr; i++)
      {
        rows_ptr[i] += rows_ptr[i - 1];
      }

    // Load the values and column numbers of transposed matrix


    for (i = 0; i < matrix.n_rows; i++)
      {
        j1 = matrix.rows_ptr[i];
        j2 = matrix.rows_ptr[i + 1];
        for (j = j1; j < j2; j++)
          {
            row_ind = matrix.cols_ind[j] - rows_offset;
            cols_ind[rows_ptr[row_ind]] = i + cols_offset;
            values[rows_ptr[row_ind]] = matrix.values[j];
            rows_ptr[row_ind]++;
          }
      }
    // rows_ptr now points to the *end* of the jth row of entries
    // instead of the beginning.  Restore rows_ptr to front of row.

    for (i = nr; i > 0; i--)
      {
        rows_ptr[i] = rows_ptr[i-1];
      }

    rows_ptr[0] = 0;

    return 0;
  }

  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::build_transpose_struct (const this_t &matrix,
      const i_type rows_offset,
      const i_type cols_offset,
      const i_type new_n_rows)
  {
    i_type i,j;
    i_type row_ind, j1, j2;
    i_type nnz;
    i_type nr;

    nnz = matrix.get_n_non_zeros ();

    if (new_n_rows == 0)
      nr = matrix.n_cols;
    else
      nr = new_n_rows;

    init (nr, matrix.n_rows, matrix.n_block_size, nnz);

    // Count the number of entries in each column of matrix (row of transposed matrix) and fill the rows_ptr array.

    for (i = 0; i < nnz; i++)
      {
        ++rows_ptr[matrix.cols_ind[i] + 1 - rows_offset];
      }

    for (i = 2; i <= n_rows; i++)
      {
        rows_ptr[i] += rows_ptr[i - 1];
      }

    // Load the values and column numbers of transposed matrix


    for (i = 0; i < matrix.n_rows; i++)
      {
        j1 = matrix.rows_ptr[i];
        j2 = matrix.rows_ptr[i + 1];
        for (j = j1; j < j2; j++)
          {
            row_ind = matrix.cols_ind[j] - rows_offset;
            cols_ind[rows_ptr[row_ind]] = i + cols_offset;
            rows_ptr[row_ind]++;
          }
      }
    // rows_ptr now points to the *end* of the jth row of entries
    // instead of the beginning.  Restore rows_ptr to front of row.

    for (i = n_rows; i > 0; i--)
      {
        rows_ptr[i] = rows_ptr[i - 1];
      }

    rows_ptr[0] = 0;

    return 0;
  }



  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::rand_init (const i_type new_n_rows,
      const i_type new_n_cols,
      const fp_type rand_value_dispersion,
      const i_type elems_in_row)
  {
    int r_code = 0;
    i_type i, j, j1, j2;
    i_type n_non_zeros;

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

  // write matrix to file with given name with or without sorting cols
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::ascii_write_in_csr_format (const std::string &file_name) const
    {
      locale_keeper lkeeper ("C", LC_NUMERIC);

      i_type j, jj, j1, j2, jj1, jj2;

      i_type b_sqr = n_block_size * n_block_size;
      FILE *fp = fopen (file_name.c_str (), "w");

      //ofstream fp(file_name.c_str ());
      if (!fp)
        //TODO: write error message
        return -1;

      fprintf (fp, "// N_ROWS\tN_COLS\tN_NON_ZEROS\tN_BLOCK_SIZE\n");
      //fp << "// N_ROWS\tN_COLS\tN_NON_ZEROS\tN_BLOCK_SIZE" << endl;

      fprintf (fp, "%d\t%d\t%d\t%d\n", n_rows, n_cols, get_n_non_zeros (), n_block_size);
      //fp << n_rows << "\t" << n_cols << "\t" << get_n_non_zeros () << "\t" << n_block_size << endl;

      fprintf (fp, "// Rows indexes[1..n_rows] (with out 0)\n");
      //fp << "// Rows indexes[1..n_rows] (with out 0)" << endl;

      for (size_t i = 1; i < rows_ptr.size (); ++i)
        {
          fprintf (fp, "%d\n", rows_ptr[i]);
          //fp << rows_ptr[i] << endl;
        }

      fprintf (fp, "// END of Rows indexes\n");
      //fp << "// END of Rows indexes" << endl;


      fprintf (fp, "// Values n_non_zeros elements\n");
      //fp << "// Values n_non_zeros elements" << endl;

      fprintf (fp, "//COLUMN\tVALUE\n");
      //fp << "//COLUMN\tVALUE" << endl;

      size_t counter = 0;
      if (rows_ptr.size ())
        {
          for (size_t i = 0; i < rows_ptr.size () - 1; ++i)
            {
              fprintf (fp, "// ROW %d\n", (int)i);
              //fp << "// ROW " << i << endl;
              j1 = rows_ptr[i];
              j2 = rows_ptr[i + 1];

              for (j = j1; j < j2; ++j, ++counter)
                {
                  fprintf (fp, "%d", cols_ind[j]);
                  //fp << cols_ind[j];// << "\t";
                  jj1 = j * b_sqr;
                  jj2 = jj1 + b_sqr;
                  for (jj = jj1; jj < jj2; ++jj)
                    fprintf (fp, "\t%.20e", values[jj]);
                    //fp << "\t" << fixed << setprecision (20) << values[jj];

                  fprintf (fp, "\n");
                  //fp<<endl;
                }
            }
        }
      if ((int)counter != get_n_non_zeros ())
        {
          BOSOUT (section::solvers, level::low) << "NO END" << "counter = " << counter << "get_n_non_zeros () = " << get_n_non_zeros () << bs_end; 
          return -1;
        }
      fprintf (fp, "// END OF FILE\n");
      //fp << "// END OF FILE" << endl;
      fflush (fp);
      fclose (fp);
      //fp.flush();
      //fp.close();
      return 0;

    }

  /*!
   * \brief read matrix from file with given name
   *
   * \param file_name -- file to read from
   *
   * \return 0 if success
   */
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::ascii_read_from_csr_format (const std::string &file_name)
  {
    locale_keeper lkeeper ("C", LC_NUMERIC);

    FILE *fp = 0;
    char buf[4096];
    int state = 0;
    i_type nc, nr, nnz = 0, nbs, b_sqr = 0;
    i_type row_ind = 0;
    i_type j_ind = 0, j = 0;
    char *start_ptr, *end_ptr;

    fp = fopen (file_name.c_str (), "r");
    BS_ASSERT (fp) (file_name);
    if (!fp)
      {
        bs_throw_exception ("Can't open file");
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

            if (init (nr, nc, nbs, nnz))
              {
                bs_throw_exception ("Can't init matrix");
                return -45;
              }

            BOSOUT << "nr = " << nr << ", nc = " << nc << ", nbs = " << nbs << ", nnz = " << nnz << bs_end;

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

    BOSOUT << "nnz = " << nnz << bs_end;

    fclose (fp);
    return 0;
  }

#ifdef _HDF5_OUT
  /*!
   * \brief dump cols_ind and values to hdf5 file
   */
  template<class fp_type, class i_type> int
  bcsr_matrix<fp_type, i_type>::dump_to_hdf5()
  {
    i_type nnz = get_n_non_zeros();

    // creating file
    H5::H5File *file = new H5::H5File(hdf5_dump_path, H5F_ACC_TRUNC);

    if (nnz != 0) // if there is data for dumping
      {
        // creating dataspace - 1-dimensional with size = nnz
        hsize_t dims[] = {nnz};
        H5::DataSpace dataspace(1, dims);

        H5::DataSet *dataset;

        // cols_ind dumping
        dataset = new H5::DataSet(file->createDataSet("cols_ind", H5::PredType::NATIVE_INT, dataspace));
        dataset->write(&cols_ind[0], H5::PredType::NATIVE_INT, dataspace);
        delete dataset;

        // values dumping
        dataset = new H5::DataSet(file->createDataSet("values", H5::PredType::NATIVE_FLOAT, dataspace));
        dataset->write(&values[0], H5::PredType::NATIVE_FLOAT, dataspace);
        delete dataset;
      }

    delete file;

    //cols_ind.clear();
    //values.clear();

    return 0;
  }

  /*!
   * \brief read hyperslabs of cols_ind and values from hdf5
   * \param hyperslab_size
   * \param cols_ind_hyperslab -- buffer for reading cols_ind
   * \param values_hyperslab -- buffer for reading values
   * \param mutex -- mutexes for synchronizing with calc thread
   * \param cond -- condition variables for synchronizing with calc thread
   * \param data_ready -- flags for synchronizing with calc thread
   */
  template<class fp_type, class i_type> int
  bcsr_matrix<fp_type, i_type>::read_from_hdf5(i_type hyperslab_size,
      i_vector_type cols_ind_hyperslab[],
      fp_vector_type values_hyperslab[],
      boost::mutex mutex[],
      boost::condition cond[],
      bool data_ready[])
  {
    i_type i, j, j1, j2, k;
    i_type nnz = get_n_non_zeros();

    H5::H5File *file = new H5::H5File(hdf5_dump_path, H5F_ACC_RDONLY);
    H5::DataSet dataset_cols_ind = file->openDataSet("/cols_ind");
    H5::DataSet dataset_values = file->openDataSet("/values");

    for (i = 0, k = 0; i < n_rows; ++i)
      {
        j1 = rows_ptr[i];
        j2 = rows_ptr[i + 1];

        for (j = j1; j < j2; j += hyperslab_size)
          {
            // preparing
            i_type real_hyperslab_size = min(hyperslab_size, j2 - j);
            hsize_t dims_memory[] = {real_hyperslab_size};
            H5::DataSpace dataspaceMemory(1, dims_memory);
            hsize_t dims_file[] = {nnz};
            H5::DataSpace dataspaceFile(1, dims_file);
            hsize_t offset[] = {j};
            hsize_t count[] = {real_hyperslab_size};
            dataspaceFile.selectHyperslab(H5S_SELECT_SET, count, offset);
            {
              // waiting while data region is blocked
              boost::mutex::scoped_lock lock(mutex[k]);
              while (data_ready[k] == true)
                cond[k].wait(lock);
              // reading data from the dump to the current data region
              cols_ind_hyperslab[k].resize(real_hyperslab_size);
              dataset_cols_ind.read(&cols_ind_hyperslab[k][0], H5::PredType::NATIVE_INT, dataspaceMemory, dataspaceFile);
              values_hyperslab[k].resize(real_hyperslab_size);
              dataset_values.read(&values_hyperslab[k][0], H5::PredType::NATIVE_FLOAT, dataspaceMemory, dataspaceFile);
              data_ready[k] = true;
            }
            cond[k].notify_one();
            // switching to the next data region
            k = (k + 1) % 2;
          }
      }

    delete file;

    return 0;
  }
#endif // #ifdef _HDF5_OUT

  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::rand_init_symm (const i_type nx, const i_type ny, const i_type nz,
      const fp_type rand_value_dispersion)
  {

    int r_code = 0;
    i_type i, j, ind[7], k;
    i_type n, nnz;
    fp_type d;

    // check input data
    if (nx < 1 || ny < 1 || nz < 1)
      // internal error
      return -1;
    n = nx * ny * nz;

    nnz = n * 7; // approx value, final nonzeros value can be equal or less

    r_code = init (n, n, 1, nnz);

    if (r_code)
      //TODO: print error message
      return -2;

    srand ((unsigned)time( NULL ));

    for (i = 0; i < nnz; ++i)
      values[i] = (fp_type)rand () / (fp_type)RAND_MAX * rand_value_dispersion;

    rows_ptr[0] = 0;
    for (i = 0, j = 0; i < n_rows; ++i)
      {
        ind[0] = i - nx * ny;
        ind[1] = i - nx;
        ind[2] = i - 1;
        ind[3] = i + 1;
        ind[4] = i + nx;
        ind[5] = i + nx * ny;
        ind[6] = i;
        d = 0;
        for (k = 0; k < 7; ++k)
          if (ind[k] >= 0 && ind[k] < n)
            {
              cols_ind[j] = ind[k];
              if (k < 6)
                {
                  values[j] = - (fp_type)rand () / (fp_type)RAND_MAX * rand_value_dispersion;
                  d += values[j];
                }
              else
                values[j] = -d + (fp_type)rand () / (fp_type)RAND_MAX * rand_value_dispersion;

              ++j;
            }
        rows_ptr[i + 1] = j;
      }
    return 0;
  }

  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::gen_2d_laplas (const i_type n)
  {
    int r_code = 0;
    i_type i, j;
    i_type N, nnz;

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

  /*!
   * \brief print staticstic of matrix
   * \return 0 if success
   */
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::print_statistic ()
  {
    //printf ("---------Matrix Statistic---------\n");
    //printf ("-- NNZ       %d\n", rows_ptr[n_rows]);
    //printf ("-- STENS     %d\n", rows_ptr[n_rows] / n_rows);
    //printf ("----------------------------------\n");
    BOSOUT (section::solvers, level::low)
      << "---------Matrix Statistic---------"
      << "-- NNZ       " << rows_ptr[n_rows]
      << "-- STENS     " << (rows_ptr[n_rows] / n_rows)
      << "----------------------------------" << bs_end;
    return 0;
  }


  /*!
   * \brief write matrix to file with given name in ij format: ROW COLUMN VALUE
   *
   * \param file_name -- file to write
   *
   * \return 0 if success
   */
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::ascii_write_in_ij_format (const std::string &file_name) const
    {
      locale_keeper lkeeper ("C", LC_NUMERIC);

      i_type i, j, j1, j2, counter;

      if (n_block_size > 1)
        {
          //printf ("N_BLOCK_SIZE > 1!\n");
          BOSERR (section::solvers, level::error) << "N_BLOCK_SIZE > 1!" << bs_end;
          return -1;
        }

      FILE *fp = fopen (file_name.c_str (), "w");
      //ofstream fp(file_name.c_str ());
      if (!fp)
        {
          //TODO: write error message
          BOSERR (section::solvers, level::error) << "Can't open file " << file_name << bs_end;
          return -1;
        }

      for (i = 0, counter = 0; i < n_rows; ++i)
        {
          j1 = rows_ptr[i];
          j2 = rows_ptr[i + 1];
          for (j = j1; j < j2; ++j, ++counter)
            {
              fprintf (fp, "%d\t%d\t%.20lf\n", i + 1, cols_ind[j] + 1, values[j]);
              //fp << (i + 1) << "\t" << (cols_ind[j] + 1) << "\t" << setprecision (20) << values[j];
            }
          fprintf (fp, "\n");
          //fp << endl;
        }
      if (counter != get_n_non_zeros ())
        {
          return -1;
        }

      fclose (fp);
      //fp.close();
      return 0;
    }


  /*!
   * check for correctness structure of matrix(rows_ptr and cols_ind)
   */
  template<class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::internal_check ()
  {
    i_type i, j, j1, j2, b_sqr, cl = 0;
    i_type i_err = 0, jj = 0, j_last = 0, r_code = 0;
    std::vector<i_type> columns;

    b_sqr = n_block_size * n_block_size;
    i_type n_non_zeros = get_n_non_zeros ();

    BOSOUT (section::solvers, level::debug) << "check_matrix: n_rows=" << n_rows << ", n_cols=" << n_cols << ", nnz=" << n_non_zeros << ", nb=" << n_block_size << bs_end;

    columns.resize (n_cols);
    memset (&columns[0], -1, sizeof (i_type) * n_cols);

    if (rows_ptr[0] != 0)
      r_code = -1;

    for (i = 0; i < n_rows && !r_code; i++)
      {
        j1 = rows_ptr[i];
        j2 = rows_ptr[i + 1];

        if (j2 <= j1)
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

        j_last = jj;
        for (j = j1; j < j2; ++j)
          {
            cl = cols_ind[j];
            if (cl < 0 || cl >= n_cols)
              {
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
      BOSOUT (section::solvers, level::debug) << "check_matrix: OK" << bs_end;
    else if (r_code == -1)
      BOSERR (section::solvers, level::error) << "Error: rows_ptr[0] must be = 0!" << bs_end;
    else if (r_code == -2)
      BOSERR (section::solvers, level::error) << boost::format ("Error in row %d: rows[i]>=rows[i+1]! (%d,%d)") % i_err % rows_ptr[i_err] % rows_ptr[i_err+1] << bs_end;
    else if (r_code == -3)
      BOSERR (section::solvers, level::debug) << "Error in row " << i_err << ": rows[i] - rows[i+1] > n_cols! (" << rows_ptr[i_err] << "," << rows_ptr[i_err+1] << ")" << bs_end;
    else if (r_code == -4)
      {
        BOSERR (section::solvers, level::debug) << "Error in row " << i_err << ": duplicate cols (" << cl << ")!" << bs_end;
      }
    else if (r_code == -5)
      BOSERR (section::solvers, level::debug) << "Error: cols_ind must be in 0 <= .. < n_cols, but" << cl << ")" << bs_end;

    fflush(stdout);
    return r_code;
  }

  /**
   * \brief Merge two BCSR matrix into one
   *
   * \param lhs Pointer to left-hand side matrix
   * \param rhs Pointer to right-hand side matrix
   * \return 0 if success
   */
  template <class fp_vector_t, class i_vector_t> int
  bcsr_matrix<fp_vector_t, i_vector_t>::merge (const this_t *lhs, const this_t *rhs, int key, const item_array_t &cfl_vector)
  {
    return merge_matrices_impl <this_t> (this, lhs, rhs, get_diag_ind (), key, cfl_vector).merge ();
  }

  ////! return copy of current matrix
  //template <class fp_vector_t, class i_vector_t> typename bcsr_matrix<fp_vector_t, i_vector_t>::sp_matrix_t
  //bcsr_matrix<fp_vector_t,  i_vector_t>::prepare_matrix () const
  //{
  //  //sp_bcsr_matrix_t sp (BS_KERNEL.create_object (this_t::bs_type ()));

  //  //sp->init(*this);
  //  //return sp;
  //
  //  return this;
  //}

  //template <typename fp_vector_t, typename i_vector_t>
  //void
  //bcsr_matrix <fp_vector_t, i_vector_t>::clear ()
  //{
  //  index_array_t tmp_rows_ptr, tmp_diag_ind, tmp_cols_ind;
  //  item_array_t tmp_values;
  //  std::swap (rows_ptr, tmp_rows_ptr);
  //  std::swap (diag_ind, tmp_diag_ind);
  //  std::swap (cols_ind, tmp_cols_ind);
  //  std::swap (values, tmp_values);
  //}
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(bcsr_matrix, (class)(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(bcsr_matrix, (class)(class));

  BLUE_SKY_TYPE_IMPL_T_EXT(2, (bcsr_matrix<seq_vector<float>, seq_vector<int> >) , 2, (matrix_base<seq_vector<float>, seq_vector<int> >), "bcsr_matrix<float, int>", "Block CSR Matrix class", "Realization of Block CSR Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(2, (bcsr_matrix<seq_vector<double>, seq_vector<int> >) , 2, (matrix_base<seq_vector<double>, seq_vector<int> >), "bcsr_matrix<double, int>", "Block CSR Matrix class", "Realization of Block CSR Matricies", false);

}//bs ns
