#ifndef PY_MATRIX_BASE_H_
#define PY_MATRIX_BASE_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include "bcsr_matrix.h"

namespace blue_sky
  {
  namespace python
    {

    // Sergey Miryanov at 07.04.2008
    // for casting into python from child classes to parent class
    // we should export base class
    template <class fp_vector_type, class i_vector_type>
    class BS_API_PLUGIN py_matrix_base : public py_objbase
      {
      public:

        typedef py_matrix_base<fp_vector_type, i_vector_type>     base_t;
        typedef matrix_base <fp_vector_type, i_vector_type>       matrix_t;
        typedef matrix_t                                          wrapped_t;
        typedef typename matrix_t::item_t                         item_t;
        typedef typename matrix_t::index_t                        index_t;
        typedef seq_vector <item_t>                               item_seq_array_t;
        typedef seq_vector <index_t>                              index_seq_array_t;

      public:

        py_matrix_base()
            : py_objbase(BS_KERNEL.create_object(matrix_t::bs_type()))
        {
        }

        py_matrix_base(sp_obj sp_obj_)
            : py_objbase(sp_obj_)
        {

        }

        py_matrix_base (const type_descriptor &td)
            : py_objbase (td)
        {

        }

        int get_n_block_size () const
        {
          return this->get_spx (this)->n_block_size;
        }
        void set_n_block_size (index_t block_size)
        {
          this->get_spx (this)->n_block_size = block_size;
        }

        void set_n_rows (index_t n_rows)
        {
          this->get_spx (this)->n_rows = n_rows;
        }
        void set_n_cols (index_t n_cols)
        {
          this->get_spx (this)->n_cols = n_cols;
        }
        index_t get_n_rows () const
        {
          return this->get_spx (this)->n_rows;
        }
        index_t get_n_cols () const
        {
          return this->get_spx (this)->n_cols;
        }

      };

  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_MATRIX_BASE_H_
