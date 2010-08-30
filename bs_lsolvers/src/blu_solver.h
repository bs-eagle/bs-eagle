#ifndef __BLU_SOLVER_H__
#define __BLU_SOLVER_H__
/*!
  \file blu_solver.h
  \brief file include declarations of functions to build LU decomposition and to solve matrix
*/
#include <string>
#include <sstream>

#include "lsolver_iface.h"
#include "lu_decomposition.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include "dens_matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  /**
   * \brief BLU solver for dens matrix
   */
  template <class strategy_t>
  class BS_API_PLUGIN blu_solver: public lsolver_iface<strategy_t>
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      //! matrix interface type
      typedef matrix_iface<strategy_t>                  matrix_t;
      //! internal fp type
      typedef typename strategy_t::fp_type_t                    fp_type_t;
      //! internal integer type
      typedef typename strategy_t::i_type_t                     i_type_t;
      //! internal matrix fp type
      typedef typename strategy_t::fp_storage_type_t            fp_storage_type_t;

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;

      //typedef linear_solver_base<strategy_t>      this_t;         ///< typedef to this type
      typedef lsolver_iface<strategy_t>           base_t;         ///< typedef to this type. in child classes used as a short name of base class

      typedef smart_ptr<matrix_t, true>           sp_matrix_t;    ///< short name to smart pointer on matrix

      typedef prop_iface<float, int, std::string, bool> prop_t;
      typedef smart_ptr<base_t, true>                   sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>               sp_prop_t;      ///< short name to smart pointer to properties holder class
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
    private:
      //! dens matrix type
      typedef dens_matrix_iface<strategy_t>  dens_matrix_iface_t;
      //! SP to dens matrix
      typedef smart_ptr<dens_matrix_iface_t, true>                                            sp_dens_matrix_t;
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      ~blu_solver ()
        {}

      virtual int solve (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol);

      virtual int solve_prec (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol);

      // setup
      virtual int setup (sp_matrix_t matrix);

      //! set preconditioner
      virtual void set_prec (sp_base_t /*prec_*/)
        {
          //prec = prec_;
        }
      
      virtual void set_prop(sp_prop_t prop_);

      //! get properties
      virtual sp_prop_t get_prop() 
        {
          return prop;
        }

      //! return final residual
      virtual fp_type_t get_final_residual () const 
        {
          return 0;
        }

      //! return number of used iterations
      virtual int get_niters () const
        {
          return 1;  // preconditioner
        }

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "Block linear solver for dense matrix using LU decomposition.\n";
          s << "Properties:\n";
          s << prop->py_str ();

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
    protected:
      void init_prop ();

    protected:


    public:
      BLUE_SKY_TYPE_DECL (blu_solver);
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
      sp_prop_t         prop;         //!< properties for solvers
      int block_size_idx;

    protected:

      blu_solver_impl<fp_type_t, i_type_t, fp_storage_type_t> impl;

    };

} // namespace blue_sky

#endif //__BLU_SOLVER_H__
