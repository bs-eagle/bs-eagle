#ifndef TWO_STAGE_PREC_H__
#define TWO_STAGE_PREC_H__
/*!
 * \file two_stage_prec.h
 * \brief class declaration for two stage preconditioner
 * \date 2006-07-26
 */
#include <string>
#include <sstream>

#ifdef _MPI
#include "mpi_csr_matrix.h"
#endif //_MPI
#include "memory_macroses.h"

#include "two_stage_prec_iface.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  template <class strategy_t>
  class BS_API_PLUGIN two_stage_prec: public two_stage_prec_iface<strategy_t>
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      //! matrix interface type
      typedef matrix_iface<strategy_t>                  matrix_t;
      //! internal fp type
      typedef typename strategy_t::fp_type_t            fp_type_t;
      typedef typename strategy_t::fp_storage_type_t    fp_storage_type_t;
      //! internal integer type
      typedef typename strategy_t::i_type_t             i_type_t;

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;

      //! prop 
      typedef prop_iface<float, int, std::string, bool>        prop_t;
      typedef lsolver_iface<strategy_t>                                 base_t;         ///< typedef to this type in child classes used as a short name of base class
      typedef smart_ptr<base_t, true>                                   sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>                                   sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>                                 sp_matrix_t;    ///< short name to smart pointer on matrix_t


      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      // destructor
      virtual ~two_stage_prec ();

      // solve
      virtual int solve (sp_matrix_t matrix, sp_fp_array_t sp_rhs, sp_fp_array_t sol);

      virtual int solve_prec (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol);

      // setup
      virtual int setup (sp_matrix_t matrix);

      //! set preconditioner
      virtual void set_prec (sp_base_t prec_)
        {
          prec = prec_;
        }

      //! set preconditioner
      virtual void set_prec_2 (sp_base_t prec_)
        {
          prec_2 = prec_;
        }
      
      virtual void set_prop (sp_prop_t prop_);

      //! get properties
      virtual sp_prop_t get_prop() 
        {
          return prop;
        }

      //! return final residual
      virtual fp_type_t get_final_residual () const 
        {
          return -1;
        }

      //! return number of used iterations
      virtual int get_niters () const
        {
          return 1;
        }

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "Two stage preconditioner.\n";
          if (prec)
            {
              s << "preconditioner 1:\n" << prec->py_str ();
            }
          if (prec_2)
            {
              s << "preconditioner 2:\n" << prec_2->py_str ();
            }
          s << "Properties:\n";
          s << prop->py_str ();

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN

    protected:
      void init_prop ();
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
    protected:
      sp_base_t         prec;         //!< pointer to the preconditioner
      sp_base_t         prec_2;         //!< pointer to the preconditioner
      sp_prop_t         prop;         //!< properties for solvers

      sp_fp_array_t     sp_r;
      sp_fp_array_t     sp_w;

      int               success_idx;

    public:
      BLUE_SKY_TYPE_DECL (two_stage_prec);
    };
} // namespace blue_sky
#endif // TWO_STAGE_PREC_H__
