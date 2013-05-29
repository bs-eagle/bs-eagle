/*!
 * \file bcsr_ilu_prec.h
 * \brief class declaration for ILU preconditioner for CSR matrix
 * \author Borschuk Oleg
 * \date 2006-11-07
 */
#ifndef BCSR_ILU__PREC__H__
#define BCSR_ILU__PREC__H__

#include <string>
#include <sstream>
#include "bs_misc.h"

#include "lsolver_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bcsr_matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  /**
   * \brief ILU preconditioner for CSR matrix
   */
  class BS_API_PLUGIN bcsr_ilu_prec: public lsolver_iface
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      //! matrix interface type
      typedef matrix_iface                                      matrix_t;
      typedef lsolver_iface                                     base_t;         ///< typedef to this type. in child classes used as a short name of base class

      typedef smart_ptr<matrix_t, true>                         sp_matrix_t;    ///< short name to smart pointer on matrix

      typedef prop_iface                                        prop_t;
      typedef smart_ptr<base_t, true>                           sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>                           sp_prop_t;      ///< short name to smart pointer to properties holder class
      //! BCSR matrix type
      typedef bcsr_matrix_iface                                 bcsr_matrix_iface_t;
      //! SP to BCSR matrix
      typedef smart_ptr<bcsr_matrix_iface_t, true>              sp_bcsr_matrix_t;
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      ~bcsr_ilu_prec ();

      virtual int solve (sp_matrix_t matrix, spv_double rhs, spv_double sol);

      virtual int solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol);

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
      virtual t_double get_final_residual () const
        {
          return 0;
        }

      //! return number of used iterations
      virtual int get_niters () const
        {
          return 1;  // preconditioner
        }

      //sp_bcsr_matrix_t get_ilu_matrix () const;

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "Incomplite LU decomposition for bloc CSR matrix.\n";
          s << "Properties:\n";
          s << wstr2str (prop->py_str ());

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
    protected:
      void init_prop ();

    protected:


    public:
      BLUE_SKY_TYPE_DECL (bcsr_ilu_prec);
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
      sp_prop_t         prop;         //!< properties for solvers
      sp_bcsr_matrix_t  lu_matrix;   //!< pointer to the LU matrix
      
    private:

    };


} // namespace blue_sky
#endif // CSR_ILU__PREC__H__
