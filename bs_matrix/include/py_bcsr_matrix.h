/**
 * @file py_bcsr_matrix.h
 * @brief
 * @author Morozov Andrey
 * @date 2008-04-02
 */
#ifndef PY_BCSR_MATRIX_H_
#define PY_BCSR_MATRIX_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include "bcsr_matrix.h"
#include "py_matrix_base.h"
#include "throw_exception.h"

namespace blue_sky
  {
  namespace python
    {

    //template <class fp_vector_type, class i_vector_type>
    //class BS_API_PLUGIN py_bcsr_matrix : public py_matrix_base<fp_vector_type, i_vector_type>
    //  {
    //  public:

    //    friend class bcsr_matrix<fp_vector_type, i_vector_type>;

    //    //types
    //    typedef typename fp_vector_type::value_type           fp_type_t;
    //    typedef typename i_vector_type::value_type            i_type_t;

    //    typedef fp_vector_type                                fp_vector_t;
    //    typedef i_vector_type                                 i_vector_t;

    //    typedef bcsr_matrix<fp_vector_type, i_vector_type>    matrix_t;
    //    typedef fp_type_t                                     item_t;
    //    typedef i_type_t                                      index_t;
    //    typedef fp_vector_t                                   item_array_t;
    //    typedef i_vector_t                                    index_array_t;

    //    typedef matrix_t                                      bcsr_matrix_t;
    //    typedef py_bcsr_matrix<fp_vector_type, i_vector_type> py_bcsr_matrix_t;
    //    typedef smart_ptr<matrix_t, true>                     sp_matrix_t;

    //    typedef py_bcsr_matrix_t                              py_this_t;
    //    typedef py_matrix_base<fp_vector_type, i_vector_type> py_matrix_base_t;
    //    typedef py_matrix_base_t                              base_t;

    //    //using base_t::template get_spx<matrix_t>;

    //    py_bcsr_matrix()
    //        : py_matrix_base_t (BS_KERNEL.create_object (matrix_t::bs_type ()))
    //    {

    //    }

    //    py_bcsr_matrix(sp_obj sp_obj_)
    //        : py_matrix_base_t (sp_obj_)
    //    {

    //    }

    //    //methods
    //    // initialize memory
    //    sp_matrix_t get_mat()
    //    {
    //      return static_cast<sp_matrix_t>(this->sp);
    //    }
    //    const sp_matrix_t get_mat() const
    //      {
    //        return static_cast<const sp_matrix_t>(this->sp);
    //      }

    //    int init (const matrix_t &matrix)
    //    {
    //      return get_mat()->init(matrix);
    //    }

    //    //// initialize memory
    //    int init (const i_type_t new_n_rows, const i_type_t new_n_cols, const i_type_t new_n_blok_size,
    //              const i_type_t new_n_non_zeros)
    //    {
    //      return get_mat()->init(new_n_rows,new_n_cols,new_n_blok_size,new_n_non_zeros);
    //    }

    //    //// initialize memory
    //    int init_struct (const i_type_t new_n_rows, const i_type_t new_n_cols, const i_type_t new_n_non_zeros)
    //    {
    //      return get_mat()->init_struct(new_n_rows,new_n_cols,new_n_non_zeros);
    //    }

    //    //// initialize indexies fo main diagonal
    //    //TODO - do smth with default params
    //    int init_diag_ind (int diag_may_not_exist = 1, int row_offset = 0)
    //    {
    //      return get_mat()->init_diag_ind(diag_may_not_exist,row_offset);
    //    }

    //    //// initialize only private arrays
    //    int alloc_rows_ptr (const i_type_t new_n_rows)
    //    {
    //      return get_mat()->alloc_rows_ptr(new_n_rows);
    //    }
    //    int alloc_cols_ind (const i_type_t new_n_non_zeros)
    //    {
    //      return get_mat()->alloc_cols_ind(new_n_non_zeros);
    //    }
    //    int alloc_values (const i_type_t new_n_non_zeros)
    //    {
    //      return get_mat()->alloc_values(new_n_non_zeros);
    //    }
    //    int alloc_cols_ind_and_values (const i_type_t new_n_non_zeros)
    //    {
    //      return get_mat()->alloc_cols_ind_and_values(new_n_non_zeros);
    //    }

    //    int rand_init (const i_type_t new_n_rows, const i_type_t new_n_cols,
    //                   const fp_type_t rand_value_dispersion, const i_type_t elems_in_row)
    //    {
    //      return get_mat()->rand_init(new_n_rows, new_n_cols, rand_value_dispersion, elems_in_row);
    //    }
    //    int rand_init_symm (const i_type_t nx, const i_type_t ny, const i_type_t nz,
    //                        const fp_type_t rand_value_dispersion)
    //    {
    //      return get_mat()->rand_init_symm(nx, ny, nz, rand_value_dispersion);
    //    }

    //    int gen_2d_laplas (const i_type_t n)
    //    {
    //      return get_mat()->gen_2d_laplas(n);
    //    }

    //    // build B = A^T from given csr matrix, return 0 if success,
    //    // rows number and offset can be set manually with new_n_rows and rows_offset
    //    int build_transpose (const matrix_t &matrix, const i_type_t rows_offset = 0,
    //                         const i_type_t cols_offset = 0, const i_type_t new_n_rows = 0)
    //    {
    //      return get_mat()->build_transpose(matrix, rows_offset, cols_offset, new_n_rows);
    //    }

    //    // build B = A^T from given csr matrix using matrix structure only, without values, return 0 if success,
    //    // rows number and offset can be set manually with new_n_rows and rows_offset
    //    int build_transpose_struct (const matrix_t &matrix,
    //                                const i_type_t rows_offset = 0,
    //                                const i_type_t cols_offset = 0,
    //                                const i_type_t new_n_rows = 0)
    //    {
    //      return get_mat()->build_transpose_struct(matrix, rows_offset, cols_offset, new_n_rows);
    //    }

    //    //////////////////////////////////////////////////////////////////////////
    //    void set_values (const item_array_t &values)
    //    {
    //      this->template get_spx <matrix_t> ()->get_values ().assign (values.begin (), values.end ());
    //    }
    //    void set_rows (const index_array_t &rows)
    //    {
    //      this->template get_spx <matrix_t> ()->get_rows_ptr ().assign (rows.begin (), rows.end ());
    //    }
    //    void set_cols (const index_array_t &cols)
    //    {
    //      this->template get_spx <matrix_t> ()->get_cols_ind ().assign (cols.begin (), cols.end ());
    //    }
    //    void set_diags (const index_array_t &diags)
    //    {
    //      this->template get_spx <matrix_t> ()->get_diag_ind ().assign (diags.begin (), diags.end ());
    //    }
    //    void set_n_rows (index_t n_rows)
    //    {
    //      this->template get_spx <matrix_t> ()->n_rows = n_rows;
    //    }
    //    void set_n_cols (index_t n_cols)
    //    {
    //      this->template get_spx <matrix_t> ()->n_cols = n_cols;
    //    }
    //    index_t get_n_rows () const
    //      {
    //        return this->template get_spx <matrix_t> ()->n_rows;
    //      }
    //    index_t get_n_cols () const
    //      {
    //        return this->template get_spx <matrix_t> ()->n_cols;
    //      }
    //    //////////////////////////////////////////////////////////////////////////

    //    //! return values array
    //    fp_vector_t &get_values ()
    //    {
    //      return get_mat()->get_values();
    //    }

    //    //! return values array for read only
    //    //TODO: should not be locked
    //    const fp_vector_t &get_values_for_read () const
    //      {
    //        return get_mat().get()->get_values ();
    //      }

    //    //! return rows_ptr array
    //    i_vector_t &get_rows_ptr ()
    //    {
    //      return get_mat()->get_rows_ptr();
    //    }

    //    //! return rows_ptr array for read only
    //    //TODO: should not be locked
    //    const i_vector_t &get_rows_ptr_for_read () const
    //      {
    //        return get_mat().get()->get_rows_ptr ();
    //      }

    //    //! return cols_ind array
    //    i_vector_t &get_cols_ind ()
    //    {
    //      return get_mat()->get_cols_ind();
    //    }

    //    //! return cols_ind array for read only
    //    //TODO: should not be locked
    //    const i_vector_t &get_cols_ind_for_read () const
    //      {
    //        return get_mat().get()->get_cols_ind ();
    //      }

    //    //! return diag_ind array
    //    i_vector_t &get_diag_ind ()
    //    {
    //      return get_mat()->get_diag_ind();
    //    }

    //    //! return diag_ind array for read only
    //    const i_vector_t &get_diag_ind_for_read () const
    //      {
    //        return get_mat().get()->get_diag_ind ();
    //      }

    //    //! return number of nonzeros elements
    //    i_type_t get_n_non_zeros () const
    //      {
    //        return get_mat()->get_n_non_zeros();
    //      }

    //    //! check for correctness the structure of matrix (rows_ptr and cols_ind)
    //    int internal_check ()
    //    {
    //      return get_mat()->internal_check();
    //    }


    //    // calculate matrix vector product, v -- input vector, r -- output vector
    //    // r += A * v
    //    // return 0 if success
    //    int matrix_vector_product (const fp_vector_t &v, fp_vector_t &r)
    //    {
    //      return get_spx <matrix_t> ()->matrix_vector_product (v, r);
    //    }

    //    // calculate matrix^t vector product, v -- input vector, r -- output vector
    //    // r += A^t * v
    //    // return 0 if success
    //    int matrix_vector_product_t (const fp_vector_t &v, fp_vector_t &r)
    //    {
    //      return get_mat()->matrix_vector_product_t(v, r);
    //    }

    //    //// calculate linear combination r = alpha * Au + beta * v
    //    //int calc_lin_comb (const fp_type_t alpha, const fp_type_t beta,
    //    //                   const fp_vector_t &u, const fp_vector_t &v, fp_vector_t &r)
    //    //{
    //    //  return get_mat()->calc_lin_comb(alpha, beta, u, v, r);
    //    //}

    //    // return total amount of allocated memory
    //    virtual i_type_t get_allocated_memory_in_bytes ()
    //    {
    //      return get_mat()->get_allocated_memory_in_bytes();
    //    }

    //    // write matrix to file with given name with or without sorting cols
    //    int ascii_write_in_csr_format (const std::string &file_name)
    //    {
    //      return get_mat()->ascii_write_in_csr_format(file_name);
    //    }

    //    // print matrix to ascii file in IJ format
    //    int ascii_write_in_ij_format (const std::string &file_name)
    //    {
    //      return get_mat()->ascii_write_in_ij_format(file_name);
    //    }

    //    // read matrix from file with given name (csr format)
    //    int ascii_read_from_csr_format (const std::string &file_name)
    //    {
    //      return get_mat()->ascii_read_from_csr_format(file_name);
    //    }

    //    int merge(py_this_t lhs, py_this_t rhs)
    //    {
    //      const sp_matrix_t &lhs_m(lhs.get_sp());
    //      const sp_matrix_t &rhs_m(rhs.get_sp());

    //      return this->template get_lspx <matrix_t> ()->merge (lhs_m, rhs_m);
    //    }

    //    int get_n_block_size () const
    //      {
    //        return this->get_spx (this)->n_block_size;
    //      }
    //    void set_n_block_size (index_t block_size)
    //    {
    //      this->template get_lspx <matrix_t> ()->n_block_size = block_size;
    //    }

    //  };

    //template <typename object_t>
    //static smart_ptr <object_t, true>
    //construct_matrix ()
    //{
    //  smart_ptr <object_t, true> obj = BS_KERNEL.create_object (object_t::bs_type ());
    //  if (!obj)
    //  {
    //    bs_throw_exception ("Can't create solver");
    //  }

    //  return obj;
    //}

    //template <typename object_t>
    //static boost::python::object 
    //init_matrix (boost::python::object py_obj)
    //{
    //  using namespace boost::python;
    //  object return_value = call_method <object> (py_obj.ptr (), "__cons__");

    //#ifdef _DEBUG
    //  object_t *obj = extract <object_t *> (py_obj);
    //  if (!obj)
    //  {
    //    bs_throw_exception ("Can't extract c++ object from python object");
    //  }
    //#endif

    //  return return_value;
    //}

    typedef bcsr_matrix <seq_vector<float>,  seq_vector<int> >    bcsr_matrix_fi;
    typedef bcsr_matrix <seq_vector<double>, seq_vector<int> >    bcsr_matrix_di;

    typedef matrix_base <seq_vector<double>, seq_vector<int> >    base_matrix_di;
    typedef matrix_base <seq_vector<float>,  seq_vector<int> >    base_matrix_fi;
  }
}//ns bs
#endif //BSPY_EXPORT_PLUGIN

#endif//PY_BCSR_MATRIX_H_
