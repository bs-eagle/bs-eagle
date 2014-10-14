/*!
 * \file two_stage_prec.cpp
 * \brief implementation of tow stage preconditioner
 * \date 2006-08-04
 */

#include "bs_common.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "two_stage_prec.h"

namespace blue_sky
  {
    //////////////////////////////////////////////////////////////////////////
    //  two_stage_prec

    //! constructor
    
    two_stage_prec::two_stage_prec (bs_type_ctor_param /*param*/)
      : two_stage_prec_iface ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
      sp_r = BS_KERNEL.create_object (v_double::bs_type ());
      sp_w = BS_KERNEL.create_object (v_double::bs_type ());
    }

    //! copy constructor
    
    two_stage_prec::two_stage_prec(const two_stage_prec &solver)
      : bs_refcounter (), two_stage_prec_iface ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    
    two_stage_prec::~two_stage_prec ()
    {}

    //! set solver's properties
    
    void two_stage_prec::set_prop (sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }

     void
    two_stage_prec::init_prop ()
      {
        prop->add_property_b (false, success_idx, L"True if solver successfully convergent");
      }
    
    int two_stage_prec::solve (sp_matrix_t matrix, 
                                        spv_double sp_rhs, 
                                        spv_double sp_sol)
    {
      BS_ERROR (matrix, "two_stage_prec_solve");
      BS_ERROR (sp_rhs->size (), "two_stage_prec_solve");
      BS_ERROR (sp_sol->size (), "two_stage_prec_solve");
      BS_ERROR (prop, "two_stage_prec_solve");

      BS_ASSERT (prec);
      BS_ASSERT (prec_2);

      t_long n = matrix->get_n_rows () * matrix->get_n_block_size ();
      BS_ASSERT ((size_t)n == sp_rhs->size ()) (n) (sp_rhs->size ());

      if (prec->solve (matrix, sp_rhs, sp_sol))
        {
          bs_throw_exception ("TWO_STAGE: PREC 1 failed");
        }

      sp_r->resize (n);
      sp_w->resize (n);
      sp_r->assign (0);
      sp_w->assign (0);

      t_double *sol = &(*sp_sol)[0];
      t_double *w = &(*sp_w)[0];

      //r_array = rhs - sol - new rhs - step 5 in Olegs book
      matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);

      // solve system with ILU - step 6
      if (prec_2->solve (matrix, sp_r, sp_w))
        {
          bs_throw_exception ("TWO_STAGE: PREC 2 failed");
        }

      // Then we make a correction of our solution x = x + w - step 7
      t_long i = 0;
      t_long n2 = n - (n % 4);
      for (; i < n2; i+=4)
        {
          sol[i]        += w[i];
          sol[i + 1]    += w[i + 1];
          sol[i + 2]    += w[i + 2];
          sol[i + 3]    += w[i + 3];
        }

      for (; i < n; i++)
        {
          sol[i] += w[i];
        }

      return 0;
    }

    
    int two_stage_prec::solve_prec (sp_matrix_t matrix, 
                                             spv_double sp_rhs, 
                                             spv_double sp_sol)
    {
      return solve (matrix, sp_rhs, sp_sol);
    }

    /**
    * @brief setup for CGS
    *
    * @param matrix -- input matrix
    *
    * @return 0 if success
    */
     int
    two_stage_prec::setup (sp_matrix_t matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("CGS: Passed matrix is null");
        }

      BS_ASSERT (prec);
      BS_ASSERT (prec_2);
      BS_ASSERT (prop);
      prec->setup (matrix);
      prec_2->setup (matrix);
      
      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (two_stage_prec);
    BLUE_SKY_TYPE_STD_COPY (two_stage_prec);

    BLUE_SKY_TYPE_IMPL (two_stage_prec, two_stage_prec_iface, "two_stage_prec", "Two stage preconditioner", "Two stage preconditioner");

} // namespace blue_sky
