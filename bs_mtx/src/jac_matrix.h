#ifndef _JAC_MATRIX_H
#define _JAC_MATRIX_H

#include <string>
#include <sstream>

#include "jac_matrix_iface.h"
#include "bdiag_matrix.h"
#include "bcsr.h"


namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  
  class BS_API_PLUGIN jac_matrix: public jac_matrix_iface
    {
    public:

      typedef bcsr_matrix_iface                                 bcsr_matrix_iface_t;
      typedef bdiag_matrix_iface                                bdiag_matrix_iface_t;

      typedef bcsr                                              bcsr_matrix_t;
      typedef bdiag_matrix                                      bdiag_matrix_t;

      typedef smart_ptr<bcsr_matrix_t, true>                    sp_bcsr_matrix_t;
      typedef smart_ptr<bdiag_matrix_t, true>                   sp_bdiag_matrix_t;

      typedef smart_ptr<bcsr_matrix_iface_t, true>              sp_bcsr_matrix_iface_t;
      typedef smart_ptr<bdiag_matrix_iface_t, true>             sp_bdiag_matrix_iface_t;

      typedef jac_matrix                                        this_t;
      typedef matrix_iface                                      base_t;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (jac_matrix);

      //! destructor
      virtual ~jac_matrix ()
      {
      }

      //-----------------------------------------
      //  matrix_iface METHODS
      //-----------------------------------------

      //! calculate matrix vector product, v -- input vector, r -- output vector
      //! r += A * v
      //! return 0 if success
      virtual int matrix_vector_product (spv_double v, spv_double r) const
      {
        return sp_accum_matrix->matrix_vector_product (v, r) 
               || sp_flux_matrix->matrix_vector_product (v, r) 
               || sp_facility_matrix->matrix_vector_product (v, r);
      }

      //! calculate matrix^t vector product, v -- input vector, r -- output vector
      //! r += A^T * v
      //! return 0 if success
      virtual int matrix_vector_product_t (spv_double v, spv_double r) const
      {
        return sp_accum_matrix->matrix_vector_product_t (v, r) 
               || sp_flux_matrix->matrix_vector_product_t (v, r) 
               || sp_facility_matrix->matrix_vector_product_t (v, r);
      }

      //! calculate linear combination r = alpha * Au + beta * v
      //! alpha, beta -- scalar
      //! v, u -- input vector
      //! r -- output vector
      //! return 0 if success
      virtual int calc_lin_comb (t_double alpha, t_double beta, spv_double u, spv_double v, spv_double r) const;

      //! return total amount of allocated memory in bytes
      virtual t_double get_allocated_memory_in_mbytes () const
      {
        return sp_accum_matrix->get_allocated_memory_in_mbytes ()
               + sp_flux_matrix->get_allocated_memory_in_mbytes ()
               + sp_facility_matrix->get_allocated_memory_in_mbytes ();
      }

      //! return block size 
      virtual t_long get_n_block_size () const
      {
        return sp_flux_matrix->get_n_block_size ();
      }

      //! return number of rows in matrix 
      virtual t_long get_n_rows () const 
      {
        return sp_flux_matrix->get_n_rows ();
      }

      //! return number of cols in matrix (return -1 in number of columns is unknown)
      virtual t_long get_n_cols () const
      {
        return sp_flux_matrix->get_n_cols ();
      }

      //! return true if number of rows is equal to the number of columns
      virtual bool is_square () const
        {
          return sp_flux_matrix->is_square ();
        }
      //! initialize vector
      virtual void init_vector (spv_double v) const
        {
          v->resize (sp_flux_matrix->get_n_rows () * sp_flux_matrix->get_n_block_size ());
          memset (&(*v)[0], 0, sizeof (t_double) * v->size ());
        }

      //-----------------------------------------
      //  jac_matrix_iface METHODS
      //-----------------------------------------

      // ------------------------------------------------
      // GET matrix private DATA methods
      // ------------------------------------------------

      //! return flux (regular) matrix
      virtual sp_bcsr_matrix_iface_t get_flux_matrix ()
        {
          return sp_flux_matrix;
        }

      //! return facility (irregular) matrix
      virtual sp_bcsr_matrix_iface_t get_facility_matrix ()
        {
          return sp_facility_matrix;
        }

      //! return accumulative matrix
      virtual sp_bdiag_matrix_iface_t get_accum_matrix ()
        {
          return sp_accum_matrix;
        }

      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      // check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const
        {
          return sp_accum_matrix->internal_check () 
                 || sp_flux_matrix->internal_check () 
                 || sp_facility_matrix->internal_check ();
        }
#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "Accumulative part of jacobian matrix:\n";
          s << sp_accum_matrix->py_str ();
          s << "Flux part of jacobian matrix:\n";
          s << sp_flux_matrix->py_str ();
          s << "Facility part of jacobian matrix:\n";
          s << sp_facility_matrix->py_str ();

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
    public:

    protected:
      
      sp_bdiag_matrix_t sp_accum_matrix;
      sp_bcsr_matrix_t sp_flux_matrix;
      sp_bcsr_matrix_t sp_facility_matrix;

    };



}//namespace blue_sky

#endif //_JACOBIAN_H

