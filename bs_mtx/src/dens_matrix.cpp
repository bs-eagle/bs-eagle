/** 
 * @file dens_matrix.cpp
 * @brief full matrix
 * @date 2009-12-08
 */
#include <sstream>

#include "dens_matrix.h" 
#include "pyublas/numpy.hpp"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  dens_matrix::dens_matrix (bs_type_ctor_param) 
        : dens_matrix_iface (),
        values (BS_KERNEL.create_object (v_float::bs_type ()))
    {
    }
  
   dens_matrix ::dens_matrix (const this_t & /*src*/) : bs_refcounter (),
        values (BS_KERNEL.create_object (v_float::bs_type ()))
     {
     }

  // TODO: 
   int
  dens_matrix::matrix_vector_product (spv_double v_, 
                                               spv_double r_) const
    {
      t_long block_rows;
      t_long block_cols;
      t_long clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      if (n_cols != (t_long)v_->size () || n_rows != (t_long)r_->size ())
        return -1;

      t_double *v = &(*v_)[0];
      t_double *r = &(*r_)[0];

      block_rows = n_rows / clb;
      if (n_rows % clb)
        ++block_rows;
      
      block_cols = n_cols / clb;
      if (n_cols % clb)
        ++block_cols;

      for (t_long i = 0; i < block_rows; ++i)
        {
          for (t_long j = 0; j < block_cols; ++j)
            {
              block_mv_product (i, j, v, r);
            }
        }
      return 0;
    }

  // TODO: 
   int
  dens_matrix::matrix_vector_product_t (spv_double v_, 
                                                 spv_double r_) const
    {
      t_long block_rows;
      t_long block_cols;
      t_long clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      if (n_cols != (t_long)r_->size () || n_rows != (t_long)v_->size ())
        return -1;

      t_double *v = &(*v_)[0];
      t_double *r = &(*r_)[0];

      block_rows = n_rows / clb;
      if (n_rows % clb)
        ++block_rows;
      
      block_cols = n_cols / clb;
      if (n_cols % clb)
        ++block_cols;

      for (t_long i = 0; i < block_rows; ++i)
        {
          for (t_long j = 0; j < block_cols; ++j)
            {
              block_mv_product_t (i, j, v, r);
            }
        }
      return 0;
    }
          
   int
  dens_matrix::calc_lin_comb (t_double alpha, 
                                       t_double beta, 
                                       spv_double u_, 
                                       spv_double v_, 
                                       spv_double r_) const
    {
      static const t_double eps = t_double (1.0e-12);
      t_long i;
      int r_code = 0;

      t_double *v = &(*v_)[0];
      t_double *r = &(*r_)[0];

      memset (r, 0, sizeof (t_double) * n_rows);

      if (fabs (alpha) > eps)
        {
          r_code = matrix_vector_product (u_, r_);
          for (i = 0; i < n_rows; ++i)
            r[i] *= alpha;
        }
      if (fabs (beta) > eps)
        {
          for (i = 0; i < n_rows; ++i)
            r[i] += v[i] * beta;
        }
      return r_code;
    }

   int
  dens_matrix::init_by_matrix (sp_dens_matrix_iface_t matrix)
    {
      return init (matrix->get_n_rows (), matrix->get_n_cols (), matrix->get_calc_block_size ());
    }

   int
  dens_matrix::init (const t_long new_n_rows, 
                              const t_long new_n_cols, 
                              const t_long block_size)
    {
      //npy_intp dims[2];
      if (new_n_rows < 1 || new_n_cols < 1)
        return -1;
      //dims[0] = n_rows = new_n_rows;
      //dims[1] = n_cols = new_n_cols;
      calc_block_size = block_size;
      values->resize (n_rows * n_cols);
      // TODO: add this method
      //values->numpy.reshape (2, dims);

      //std::fill (values->numpy.as_ublas ().data ().begin (), values->numpy.as_ublas ().data ().end (), 0);
      std::fill (values->begin (), values->end (), 0);
      return 0;
    }

   int
  dens_matrix::copy (sp_dens_matrix_iface_t matrix)
    {
      if (init (matrix->get_n_rows (), matrix->get_n_cols (), matrix->get_calc_block_size ()))
        return -3;
      values = matrix->get_values ()->clone ();
      //values->numpy.resize (matrix->get_values ()->numpy.size1 (), matrix->get_values ()->numpy.size2 ());
      //values->numpy.as_ublas () = matrix->get_values ()->numpy.as_ublas (); 
                              //                  matrix->get_values ()->numpy.as_ublas ().data ().end ());
       

      //memcpy (&values[0], &(matrix->get_values_const ())[0], n_rows * n_cols * sizeof (t_float));
      return 0;
    }

   void
  dens_matrix::block_mv_product (t_long row_block, 
                                          t_long col_block, 
                                          const t_double *v, 
                                          t_double *r) const
    {
      t_long block_n_rows;
      t_long block_n_cols;
      const t_float *block;
      const t_double *block_v;
      t_double *block_r;
      t_long clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      // calculate block start position
      //block = &values[row_block * clb * n_cols + col_block * clb];
      block = &(*values)[0] + row_block * clb * n_cols + col_block * clb;
      block_v = v + col_block * clb;
      block_r = r + row_block * clb;
      
      if ((row_block + 1) * clb <= n_rows)
        block_n_rows = clb;
      else
        block_n_rows = n_rows - row_block * clb;

      if ((col_block + 1) * clb <= n_cols)
        block_n_cols = clb;
      else
        block_n_cols = n_cols - col_block * clb;

      for (t_long i = 0; i < block_n_rows; ++i)
        {
          for (t_long j = 0; j < block_n_cols; ++j)
            {
              block_r[i] += block[i * n_rows + j] * block_v[j];
            }
        }
    }
   void
  dens_matrix::block_mv_product_t (t_long row_block, 
                                            t_long col_block, 
                                            const t_double *v, 
                                            t_double *r) const
    {
      t_long block_n_rows;
      t_long block_n_cols;
      const t_float *block;
      const t_double *block_v;
      t_double *block_r;
      t_long clb = calc_block_size < 1 ? (n_rows + n_cols) : calc_block_size; 

      // calculate block start position
      block = &(*values)[0] + row_block * clb * n_cols + col_block * clb;
      block_v = v + row_block * clb;
      block_r = r + col_block * clb;
      
      if ((row_block + 1) * clb <= n_rows)
        block_n_rows = clb;
      else
        block_n_rows = n_rows - row_block * clb;

      if ((col_block + 1) * clb <= n_cols)
        block_n_cols = clb;
      else
        block_n_cols = n_cols - col_block * clb;

      for (t_long i = 0; i < block_n_rows; ++i)
        {
          for (t_long j = 0; j < block_n_cols; ++j)
            {
              block_r[j] += block[i * n_rows + j] * block_v[i];
            }
        }
    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string
  dens_matrix::py_str () const
    {
      stringstream s;

      s << "Dens matrix\n--------------------------------\n";
      s << "\tRows:                  " << n_rows << "\n";
      s << "\tColumns:               " << n_cols << "\n";
      s << "\tCalc Block:            " << calc_block_size << "\n";
      s << "\tAllocated memory (Mb): " << get_allocated_memory_in_mbytes () << "\n";
      s << "--------------------------------\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (dens_matrix);
  BLUE_SKY_TYPE_STD_COPY (dens_matrix);

  BLUE_SKY_TYPE_IMPL (dens_matrix, dens_matrix_iface,  "dens_matrix", "Dens Matrix class", "Realization of Dens Matricies");
}  // blue_sky namespace
