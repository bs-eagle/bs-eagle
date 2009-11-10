/**
 * \file read_b_matrix.h
 * \brief read b_matrix from file and convert this matrix to bcsr_matrix
 * \author Sergey Miryanov
 * \date 07.06.2008
 * */
#ifndef BS_READ_B_MATRIX_H_
#define BS_READ_B_MATRIX_H_

#include "bcsr_matrix.h"

#include "strategies.h"

#include "bs_kernel.h"
#include "naive_file_reader.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include "py_bcsr_matrix.h"
#endif


namespace blue_sky
  {
  namespace tools
    {

    template <class matrix_t>
    static smart_ptr <matrix_t, true >
    read_b_matrix_from_file (const char *filename, bool summ_diag)
    {
      BS_ASSERT (filename);

      typedef smart_ptr <matrix_t, true> sp_matrix_t;
      typedef typename matrix_t::item_array_t item_array_t;
      typedef typename matrix_t::index_array_t index_array_t;
      typedef typename item_array_t::value_type   item_t;
      typedef typename index_array_t::value_type  index_t;

      sp_matrix_t matrix = BS_KERNEL.create_object (matrix_t::bs_type ());
      BS_ERROR (matrix, "read_b_matrix_from_file");
      const sp_matrix_t &locked_matrix (matrix);

      index_t n_rows = 0, n_of_diag_band = 0, n_block_size = 0;
      item_array_t main_diagonal;
      item_array_t main_diagonal_acc;
      index_array_t bands_ind;
      item_array_t bands;

      naive_file_reader (filename)
      .locate_section ("N_ROWS")
      .read_item (n_rows).read_item (n_of_diag_band).read_item (n_block_size)
      .locate_section ("MAIN DIAGONAL")
      .read_list (main_diagonal)
      .locate_section ("MAIN DIAGONAL ACCUMULATIVE")
      .read_list (main_diagonal_acc)
      .locate_section ("BANDS INDEX")
      .read_list (bands_ind)
      .locate_section ("BANDS")
      .read_list (bands)
      ;

      index_t b_sqr = n_block_size * n_block_size;
      index_t n = n_rows * b_sqr;

      BS_ASSERT (n_rows && n_of_diag_band && n_block_size)(n_rows)(n_of_diag_band)(n_block_size);
      BS_ASSERT (main_diagonal.size () == (size_t)n) (main_diagonal.size ()) (n);
      BS_ASSERT (main_diagonal_acc.size () == (size_t)n) (main_diagonal_acc.size ()) (n);

      // fill bcsr_matrix
      index_array_t &rows = locked_matrix->get_rows_ptr ();
      rows.assign (n_rows + 1, 0);

      index_t value_count = 0;
      for (index_t i = 0; i < n_rows; ++i)
        {
          rows[i] += value_count;
          ++value_count;
          for (index_t j = 0, jcnt = 2 * n_of_diag_band; j < jcnt; ++j)
            {
              if (bands_ind[i * jcnt + j] >= 0)
                {
                  value_count ++;
                }
            }
        }
      rows[n_rows]  = value_count;

      index_array_t &cols   = locked_matrix->get_cols_ind ();
      item_array_t &values  = locked_matrix->get_values ();
      index_array_t &diags   = locked_matrix->get_diag_ind ();

      assign (diags, n_rows, -1);
      assign (cols, value_count, -1);
      assign (values, value_count * b_sqr, 0);

#ifdef _DEBUG
      item_t *values_ptr_ = &values[0];
      index_t *cols_ptr_  = &cols[0];
      index_t *rows_ptr_  = &rows[0];
#endif

      if (summ_diag)
        {
          // sum elements of main_diagonal and accumulative main diagonal
          for (index_t i = 0; i < n_rows; ++i)
            {
              for (index_t k = 0; k < b_sqr; ++k)
                {
                  main_diagonal [i * b_sqr + k] += main_diagonal_acc [i * b_sqr + k];
                }
            }
        }

      int marker = 0;
      for (index_t i = 0; i < n_rows; ++i)
        {
#ifndef _DEBUG
          memcpy (&values[marker * b_sqr], &main_diagonal[i * b_sqr], sizeof (item_t) * b_sqr);
#else
          for (int k = 0, f = 0, b_sqr_ = b_sqr; k < b_sqr; ++k, ++f)
            {
              BS_ASSERT (values [marker * b_sqr_ + f] == 0);
              values [marker * b_sqr_ + f] = main_diagonal [i * b_sqr + k];
            }
#endif
          BS_ASSERT (cols[marker] == -1) (cols[marker]);
          cols[marker] = i;
          diags[i] = marker;
          ++marker;

          for (index_t j = 0, jcnt = 2 * n_of_diag_band; j < jcnt; ++j)
            {
              int l = i * jcnt + j;
              index_t col = bands_ind[l];
              if (col >= 0)
                {
                  BS_ASSERT (col < n_rows)(col)(n_rows);
                  BS_ASSERT (col != i)(col)(i);

#ifndef _DEBUG
                  item_t *dst = &values[marker * b_sqr];
                  item_t *src = &bands[l * b_sqr];
                  memcpy (dst, src, sizeof (item_t) * b_sqr);
#else
                  for (int k = 0, f = 0, b_sqr_ = b_sqr; k < b_sqr; ++k, ++f)
                    {
                      BS_ASSERT (values [marker * b_sqr_ + f] == 0);
                      values [marker * b_sqr_ + f] = bands [l * b_sqr + k];
                    }
#endif
                  BS_ASSERT (cols[marker] == -1) (cols[marker]);
                  cols[marker] = col;
                  ++marker;
                }
            }
        }

      BS_ASSERT ((size_t)marker == cols.size ())(marker)(cols.size ());

      locked_matrix->n_rows       = n_rows;//(index_t)rows.size () - 1;
      locked_matrix->n_block_size = n_block_size;
      locked_matrix->n_cols       = n_rows;//(index_t)cols.size ();

      return matrix;
    }
  } // namespace tools

#ifdef BSPY_EXPORTING_PLUGIN

  namespace python {

    inline void
    py_export_read_b_matrix ()
    {
      using namespace boost::python;

      def ("read_b_matrix_fi", tools::read_b_matrix_from_file <base_strategy_fi::csr_matrix_t>);
      def ("read_b_matrix_di", tools::read_b_matrix_from_file <base_strategy_di::csr_matrix_t>);
    }

  } // namespace python
#endif

} // namespace blue_sky


#endif  // #ifndef BS_READ_B_MATRIX_H_
