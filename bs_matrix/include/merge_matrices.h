/**
 * \file merge_matrices.h
 * \brief merge two matrices into one
 * \author Sergey Miryanov
 * \date 13.01.2009
 * */
#ifndef BS_MERGE_MATRICES_H_
#define BS_MERGE_MATRICES_H_

#include "matrix_sum_spec.h"

#define CFL_NUM 1

namespace blue_sky {

  template <typename matrix_t>
  struct merge_matrices_impl
  {
    typedef typename matrix_t::index_t          index_t;
    typedef typename matrix_t::item_t           item_t;
    typedef typename matrix_t::item_array_t     item_array_t;
    typedef typename matrix_t::index_array_t    index_array_t;


    merge_matrices_impl (matrix_t *dst, const matrix_t *reg, const matrix_t *irreg, index_array_t &markers, int key_, const item_array_t &cfl_)
    : dst (dst),
    rows (dst->get_rows_ptr ()),
    cols (dst->get_cols_ind ()),
    values (dst->get_values ()),
    reg_rows (reg->get_rows_ptr ()),
    reg_cols (reg->get_cols_ind ()),
    reg_values (reg->get_values ()),
    irreg_rows (irreg->get_rows_ptr ()),
    irreg_cols (irreg->get_cols_ind ()),
    irreg_values (irreg->get_values ()),
    markers (markers),
    cfl (cfl_),
    marker (0),
    n ((index_t)rows.size ()),
    nb (reg->n_block_size),
    b_sqr (nb * nb),
    block_size_ (sizeof (item_t) * b_sqr),
    counter (0),
    key (key_)
    {
      dst->n_block_size = nb;
    }

    index_t
    merge (bool merge_vals = 1)
    {
      if (reg_rows.size () != irreg_rows.size ())
        return -1;

      bool ok = false;
      if (!reg_rows.empty () && !irreg_rows.empty ())
        {
          n  = (index_t)reg_rows.size ();
          ok = !dst->alloc_rows_ptr ((index_t) (reg_rows.size () - 1), -1)
            && lhs_rhs_stage1 ()
            && (merge_vals ? !dst->alloc_cols_ind_and_values (rows.back ()) : !dst->alloc_cols_ind (rows.back ()))
            && set_markers ()
            && stage2 (merge_vals);
        }
      else if (!reg_rows.empty ())
        {
          n  = (index_t)reg_rows.size ();
          ok = !dst->alloc_rows_ptr ((index_t) (reg_rows.size () - 1), -1)
            && lhs_stage1 ()
            && (merge_vals ? !dst->alloc_cols_ind_and_values (rows.back ()) : !dst->alloc_cols_ind (rows.back ()))
            && set_markers ()
            && stage2 (merge_vals);
        }
      else if (!irreg_rows.empty ())
        {
          n  = (index_t)irreg_rows.size ();
          ok = !dst->alloc_rows_ptr ((index_t) (irreg_rows.size () - 1), -1)
            && rhs_stage1 ()
            && (merge_vals ? !dst->alloc_cols_ind_and_values (rows.back ()) : !dst->alloc_cols_ind (rows.back ()))
            && set_markers ()
            && stage2 (merge_vals);
        }

      dst->n_rows = (index_t)rows.size () - 1;
      dst->n_cols = dst->n_rows;

      return ok ? 0 : -1;
    }

    index_t
    merge_stage1_init_rows ()
    {
      if (reg_rows.size () != irreg_rows.size ())
        return -1;

      bool ok = false;
      if (!reg_rows.empty () && !irreg_rows.empty ())
        {
          n  = reg_rows.size ();
          ok = lhs_rhs_stage1 ();
        }
      else if (!reg_rows.empty ())
        {
          n  = reg_rows.size ();
          ok = lhs_stage1 ();
        }
      else if (!irreg_rows.empty ())
        {
          n  = irreg_rows.size ();
          ok = rhs_stage1 ();
        }

      return ok ? 0 : -1;
    }

  private:

    bool
    stage2 (bool merge_vals = 1)
    {
      if (!reg_rows.empty () && !irreg_rows.empty ())
        {
          if (nb == 1)
            lhs_rhs_stage2 <sum_impl <1> > (merge_vals);
          else if (nb == 2)
            lhs_rhs_stage2 <sum_impl <2> > (merge_vals);
          else if (nb == 3)
            lhs_rhs_stage2 <sum_impl <3> > (merge_vals);
        }
      else if (!reg_rows.empty ())
        {
          if (nb == 1)
            lhs_stage2 <sum_impl <1> > (merge_vals);
          else if (nb == 2)
            lhs_stage2 <sum_impl <2> > (merge_vals);
          else if (nb == 3)
            lhs_stage2 <sum_impl <3> > (merge_vals);
        }
      else if (!irreg_rows.empty ())
        {
          if (nb == 1)
            rhs_stage2 <sum_impl <1> > (merge_vals);
          else if (nb == 2)
            rhs_stage2 <sum_impl <2> > (merge_vals);
          else if (nb == 3)
            rhs_stage2 <sum_impl <3> > (merge_vals);
        }

      return true;
    }

    bool
    lhs_rhs_stage1 ()
    {
      for (index_t i = 0; i < n - 1; i++)
        {

          index_t count = init_rows (i, reg_rows, reg_cols, irreg_rows);
          count += init_rows (i, irreg_rows, irreg_cols, irreg_rows);
          // uchet diagonalnyih elementov
          if (key == 0) counter += ((count != 0) ? 1 : 0);
          rows[i + 1] = counter;
        }

      return true;
    }
    bool
    lhs_stage1 ()
    {
      for (index_t i = 0; i < n - 1; i++)
        {
          index_t count = init_rows (i, reg_rows, reg_cols, irreg_rows);
          if (key == 0) counter += ((count != 0) ? 1 : 0);
          rows[i + 1] = counter;
        }

      return true;
    }
    bool
    rhs_stage1 ()
    {
      for (index_t i = 0; i < n - 1; i++)
        {
          index_t count = init_rows (i, irreg_rows, irreg_cols, irreg_rows);
          if (key == 0) counter += ((count != 0) ? 1 : 0);
          rows[i + 1] = counter;
        }

      return true;
    }

    index_t
    init_rows (index_t i, const index_array_t &rows_, const index_array_t &cols_, const index_array_t &irreg_rows_)
    {
      index_t j = rows_[i];
      index_t j2 = rows_[i + 1];
      index_t count = j2 - j;
      index_t cl;
      index_t cl_1 = irreg_rows_[i + 1] - irreg_rows_[i];

      static int COUNT = 0;

       if (key == 0)
         {
            for (; j < j2; j++)
              {
                index_t col = cols_[j];
                 if (markers[col] < i)
                  {
                    if (col != i)
                      {
                        counter++;
                      }
                    markers[col] = i;
                  }
             }
        }

        if (key == 1)
          {

            for (; j < j2; j++)
              {
                index_t col = cols_[j];
                BS_ASSERT (col >= 0) (col) (j) (i);
                if (markers[col] < i)
                  {
                    //if (col != i)
                      {
                        cl = irreg_rows[col + 1] - irreg_rows[col];
                        bool do_ = false;
                        if ((cl > 0 || fabs(cfl[col]) > CFL_NUM) && (cl_1 > 0 || fabs(cfl[i]) > CFL_NUM))
                          {
                            do_ = true;
                          }
                        if (col == i && fabs(cfl[i]) <= CFL_NUM && cl < 1)
                          {
                            do_ = true;
                          }

                        if (do_) counter++;

                      }
                    markers[col] = i;

                  }

              }

          }

      COUNT = 0;
      return count;
    }

    template <typename sum_impl_t>
    void
    lhs_rhs_stage2 (bool merge_vals = 1)
    {
      for (index_t i = 0; i < n - 1; i++)
        {
          index_t lhs_j		= reg_rows[i];
          index_t lhs_j2	= reg_rows[i + 1];
          index_t rhs_j		= irreg_rows[i];
          index_t rhs_j2	= irreg_rows[i + 1];
          index_t cl;
          index_t cl_1 = irreg_rows[i + 1] - irreg_rows[i];
          index_t col_1, col_2;

          //set_diag (i);

          if (key == 0)
          {

            for (; lhs_j < lhs_j2; lhs_j++)
              {
                col_2 = reg_cols[lhs_j] + 1;
                col_1 = reg_cols[lhs_j];
                cl = irreg_rows[col_2] - irreg_rows[col_1];

                stage2_inner_cycle <sum_impl_t> (i, lhs_j, reg_values, reg_cols, merge_vals);
              }

            for (; rhs_j < rhs_j2; rhs_j++)
              {
                col_2 = irreg_cols[rhs_j] + 1;
                col_1 = irreg_cols[rhs_j];
                cl = irreg_rows[col_2] - irreg_rows[col_1];

                stage2_inner_cycle <sum_impl_t> (i, rhs_j, irreg_values, irreg_cols, merge_vals);
              }
          }

          if (key == 1)
            {
              for (; lhs_j < lhs_j2; lhs_j++)
                {
                  col_2 = reg_cols[lhs_j] + 1;
                  col_1 = reg_cols[lhs_j];
                  cl = irreg_rows[col_2] - irreg_rows[col_1];
                  bool do_ = false;
                  if ((cl > 0 || fabs(cfl[col_1]) > CFL_NUM) && (cl_1 > 0 || fabs(cfl[i]) > CFL_NUM))
                    {
                      do_ = true;
                    }
                  if (i == col_1 && fabs(cfl[i]) <= CFL_NUM && cl < 1)
                    do_ = true;

                  if (do_)
                    stage2_inner_cycle <sum_impl_t> (i, lhs_j, reg_values, reg_cols, merge_vals);


                }

              for (; rhs_j < rhs_j2; rhs_j++)
                {
                  col_2 = irreg_cols[rhs_j] + 1;
                  col_1 = irreg_cols[rhs_j];
                  cl = irreg_rows[col_2] - irreg_rows[col_1];

                  bool do_ = false;
                  if ((cl > 0 || fabs(cfl[col_1]) > CFL_NUM) && (cl_1 > 0 || fabs(cfl[i]) > CFL_NUM))
                    {
                      do_ = true;
                    }
                  if (i == col_1 && fabs(cfl[i]) <= CFL_NUM && cl < 1)
                    do_ = true;

                  if (do_)
                    stage2_inner_cycle <sum_impl_t> (i, rhs_j, irreg_values, irreg_cols, merge_vals);

                }
            }
        }
    }
    template <typename sum_impl_t>
    void
    lhs_stage2 (bool merge_vals = 1)
    {
      for (index_t i = 0; i < n - 1; i++)
        {
          index_t lhs_j		= reg_rows[i];
          index_t lhs_j2	= reg_rows[i + 1];

          set_diag (i);
          for (; lhs_j < lhs_j2; lhs_j++)
            {
              stage2_inner_cycle <sum_impl_t> (i, lhs_j, reg_values, reg_cols, merge_vals);
            }
        }
    }
    template <typename sum_impl_t>
    void
    rhs_stage2 (bool merge_vals = 1)
    {
      for (index_t i = 0; i < n - 1; i++)
        {
          index_t rhs_j		= irreg_rows[i];
          index_t rhs_j2	= irreg_rows[i + 1];

          set_diag (i);
          for (; rhs_j < rhs_j2; rhs_j++)
            {
              stage2_inner_cycle <sum_impl_t> (i, rhs_j, irreg_values, irreg_cols, merge_vals);
            }
        }
    }

    void
    set_diag (index_t i)
    {
#ifdef _DEBUG
      for (index_t j = 0, cnt = marker * b_sqr; j < b_sqr; ++j)
        {
          BS_ASSERT (j + cnt < (index_t)values.size ()) (j) (cnt) (values.size ()) (marker) (b_sqr);
          values[j + cnt] = 0;
        }
#else
      memset (&values[marker * b_sqr], 0, block_size_);
#endif
      markers[i] = marker;
      ++marker;
    }

    bool
    set_markers ()
    {
      assign (markers, -1);
      marker = 0;

      return true;
    }

    /**
     * \brief Inner cycle of merging two matrices
     *
     * \param src Array of values from source matrix
     * \param src_cols Array of source column indexes
     */
    template <typename sum_impl_t>
    void
    stage2_inner_cycle (index_t i, index_t j, const item_array_t &values_, const index_array_t &cols_, bool merge_vals = 1)
    {
      BS_ASSERT (cols_.size () > (size_t)j) (j) (cols_.size ());
      BS_ASSERT (markers.size () > (size_t)(cols_[j])) (cols_[j]);
      index_t col = cols_[j];

      if (markers[col] >= rows[i])
        {
          if (merge_vals)
            {
              index_t index_dst = markers[col] * b_sqr;
              index_t index_src = j * b_sqr;

              BS_ASSERT (values.size () > (size_t)index_dst) (index_dst);
              BS_ASSERT (values_.size () > (size_t)index_src) (index_src);
              if (values.size () <= (size_t)index_dst)
                {
                  throw bs_exception ("stage2_inner_cycle", "index_dst out of range");
                }
              if (values_.size () <= (size_t)index_src)
                {
                  throw bs_exception ("stage2_inner_cycle", "index_src out of range");
                }

              item_t *pdst        = &values[index_dst];
              const item_t *psrc  = &values_[index_src];

              sum_impl_t::sum (psrc, pdst);
            }
          cols[markers[col]] = col;
        }
      else
        {
          if (merge_vals)
            {
              index_t index_dst = marker * b_sqr;
              index_t index_src = j * b_sqr;
              BS_ASSERT (values.size () > (size_t)index_dst) (index_dst);
              BS_ASSERT (values_.size () > (size_t)index_src) (index_src);
              if (values.size () <= (size_t)index_dst)
                {
                  throw bs_exception ("stage2_inner_cycle", "index_dst out of range");
                }
              if (values_.size () <= (size_t)index_src)
                {
                  throw bs_exception ("stage2_inner_cycle", "index_src out of range");
                }

              item_t *pdst        = &values[index_dst];
              const item_t *psrc  = &values_[index_src];

              sum_impl_t::copy (psrc, pdst);
            }
          markers[col] = marker;
          cols[marker] = col;
          ++marker;
        }
    }

    matrix_t              *dst;
    index_array_t         &rows;
    index_array_t         &cols;
    item_array_t          &values;
    const index_array_t   &reg_rows;
    const index_array_t   &reg_cols;
    const item_array_t    &reg_values;
    const index_array_t   &irreg_rows;
    const index_array_t   &irreg_cols;
    const item_array_t    &irreg_values;
    index_array_t         &markers;
    const item_array_t    &cfl;

    index_t               marker;
    index_t               n;
    index_t               nb;
    index_t               b_sqr;
    index_t               block_size_;
    index_t               counter;
    index_t               key;
  };
}


#endif  // #ifndef BS_MERGE_MATRICES_H_
