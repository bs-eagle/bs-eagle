#ifndef LIN_SOLVER__H
#define LIN_SOLVER__H
/*!
* \file linear_solvers.h
* \brief linear solvers declaration
* \author Borschuk Oleg
* \date 2006-07-26
*/

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bcsr_matrix.h"
#include "named_pbase_access.h"
#include "bos_report.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  /**
  * \brief properties for linear solvers
  */
  class BS_API_PLUGIN linear_solver_prop : public named_pbase
    {
      typedef property_base base_t;                 ///< short name for base class

    public:
      typedef double fp_type_t;                        ///< short name for used floating point type

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~linear_solver_prop ();

      //! set maximum number of iterations
      int set_max_iters (int n_iters);

      //! set tolerance
      void set_tolerance (fp_type_t new_tol)
      {
        if (new_tol > 10e-16)
          set_param(FP_TOLERANCE,(fp_type_t)new_tol);
      }

      //! set tolerance
      void set_matbal_tolerance (fp_type_t new_tol)
      {
        if (new_tol > 10e-16)
          set_param(FP_MATBAL_TOLERANCE,(fp_type_t)new_tol);
      }

      //! set number of iters
      void set_iters (int n_iters)
      {
        set_param(I_ITERS,n_iters);
      }

      //! set successively converged flag
      void set_success (int success)
      {
        set_param(I_SUCCESS,success);
      }

      //! set final resid
      void set_final_resid (fp_type_t final_resid)
      {
        set_param(FP_FINAL_RESID,(fp_type_t)final_resid);
      }

      //! set relative factor
      void set_relative_factor (fp_type_t relative_factor)
      {
        set_param(FP_RELATIVE_FACTOR,(fp_type_t)relative_factor);
      }

      //! return != 0 if method successfully converged
      int check_convergence () const
        {
          return get_int(I_SUCCESS);
        }

      //! return number of iteration
      int get_iters () const
        {
          return get_int(I_ITERS);
        }

      //! return relative residual denominator
      fp_type_t get_relative_factor () const
        {
          return get_float(FP_RELATIVE_FACTOR);
        }

      //! return maximum allowed number of iterations
      int get_max_iters () const
        {
          return get_int(I_MAX_ITERS);
        }

      //! return tolerance
      fp_type_t get_tolerance () const
        {
          return get_float(FP_TOLERANCE);
        }

      //! get successively converged flag
      int get_success () const
        {
          return get_int(I_SUCCESS);
        }

      //! get final resid
      fp_type_t get_final_resid () const
        {
          return get_float(FP_FINAL_RESID);
        }

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
      PROP_BASE_IDX_DECL_BEGIN(linear_solver_prop,property_base)
      FP_RELATIVE_FACTOR,               //!< Parameter index: denominator for resid array to get relative residual
      FP_TOLERANCE,                     //!< Parameter index: tolerance
      FP_MATBAL_TOLERANCE,
      FP_FINAL_RESID,                   //
      FP_TOTAL,
      I_ITERS,                          //!< Parameter index: number of iteration spend for convergence
      I_MAX_ITERS,                      //!< Parameter index: maximum number of iteration
      I_SUCCESS,                        //!< Parameter index: != 0 if successively converged
      I_TOTAL,
      LS_TOTAL,
      PROP_BASE_IDX_DECL_END

    public:
      BLUE_SKY_TYPE_DECL (linear_solver_prop);
      PBASE_ACCESS_MS(linear_solver_prop)

    public:
      virtual const std::string &get_params_name (idx_type idx);
    };

  /**
  * \brief base interface class for linear solvers
  */
  template <class strategy_t>
  class BS_API_PLUGIN linear_solver_base : public bs_node
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      typedef typename strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
      typedef typename strategy_t::item_t         item_t;         ///< short name to array item type (old array_item_t)
      typedef typename strategy_t::index_t        index_t;        ///< short name to matrix's index type
      typedef typename strategy_t::item_array_t   item_array_t;   ///< short name to array of items type (old array_t)
      typedef typename strategy_t::index_array_t  index_array_t;  ///< short name to array type
      typedef typename strategy_t::wksp_t         wksp_t;

      typedef typename strategy_t::rhs_item_t	      rhs_item_t;         ///< short name for type of rhs
      typedef typename strategy_t::rhs_item_array_t	rhs_item_array_t;   ///< short name for rhs array type

      typedef linear_solver_base<strategy_t>      this_t;         ///< typedef to this type
      typedef linear_solver_base<strategy_t>      base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<this_t, true>             sp_this_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<linear_solver_prop, true> sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>           sp_matrix_t;    ///< short name for smart pointer on matrix

      typedef strategy_t                          strategy_type;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      // constructor
      linear_solver_base (bs_type_ctor_param param, bs_node::sp_node node);
      // destructor
      virtual ~linear_solver_base ();

      //! solve, used vectors to array
      virtual int solve (matrix_t * /*matrix*/, typename strategy_t::rhs_item_array_t & /*rhs*/, typename strategy_t::item_array_t & /*sol*/)
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

      //! solve, used vectors to array
      virtual int solve_prec (matrix_t * /*matrix*/, typename strategy_t::item_array_t & /*rhs*/, typename strategy_t::item_array_t & /*sol*/)
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

      //! setup
      virtual int setup (matrix_t * /*matrix*/)
      {
        BS_ASSERT (false && "PURE CALL");
        return 0;
      }

      //! set preconditioner
      void set_prec (const this_t *new_prec)
      {
        prec = new_prec;
        set_subnode_in_tree ("prec", sp_this_t (new_prec));
      }

      //! set solver's properties
      void set_prop(const linear_solver_prop *new_prop)
      {
        prop = new_prop;
        set_subnode_in_tree ("prop", sp_prop_t (new_prop));
      }

      //! get properties
      sp_prop_t get_prop() const
        {
          return prop;
        }

    protected:

      //! set subnode in tree by name
      void set_subnode_in_tree (const std::string &name, const sp_obj &obj)
      {
        //bs_node::erase (name);
        //bs_node::insert (obj, name, false);
      }

      /**
       * \brief Trait class for control of adding new nodes to solver
       */
      struct solver_trait : bs_node::sort_traits
        {
          /**
           * \brief Class for sorting sub nodes of solver node
           */
          struct solver_key : bs_node::sort_traits::key_type
            {
              virtual bool sort_order (const key_ptr & ) const
                {
                  return true;
                }
            };

          virtual const char * sort_name () const
            {
              return "solver trait";
            };

          virtual key_ptr key_generator (const sp_link& /*l*/) const
            {
              return new solver_key ();
            }

          /**
           * \brief Check passed parameter on accordance to allowed object types
           *
           * \param l Link on new object
           * \return True if parameter allowed, false in other cases
           */
          virtual bool accepts (const sp_link& l)
          {
            if (l->name() == "prec")
              {
                sp_this_t sp(l->data (), bs_dynamic_cast ());
                return sp;
              }
            else if (l->name() == "prop")
              {
                sp_prop_t sp (l->data (), bs_dynamic_cast ());
                return sp;
              }

            return false;
          }
        };

    private:
      void init();

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    private:


    protected:
      wksp_t            wksp;         //!< workspace array
      sp_this_t         prec;         //!< pointer to the preconditioner
      sp_prop_t         prop;         //!< properties for solvers

    public:
      BLUE_SKY_TYPE_DECL (linear_solver_base<strategy_t>);
    };

  /**
  * @brief GMRES linear solver
  */
  template <class strategy_t>
  class BS_API_PLUGIN gmres_solver2 : public linear_solver_base<strategy_t>
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      typedef typename strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
      typedef typename strategy_t::item_array_t   item_array_t;   ///< short name to array type
      typedef typename strategy_t::item_t         item_t;         ///< short name to array item type
      typedef typename strategy_t::index_t        index_t;        ///< short name to matrix's index type
      typedef typename strategy_t::index_array_t  index_array_t;
      typedef typename strategy_t::barrier_t      barrier_t;

      typedef typename strategy_t::rhs_item_t	  rhs_item_t;     ///< short name for type of rhs
      typedef typename strategy_t::rhs_item_array_t	rhs_item_array_t;   ///< short name for rhs array type

      typedef bcsr_matrix<rhs_item_array_t, index_array_t>    bcsr_matrix_t;        ///< short name for used matrix
      typedef smart_ptr<bcsr_matrix_t, true>                  sp_bcsr_matrix_t;     ///< short name for smart_pointer on used matrix

      typedef linear_solver_base<strategy_t>      this_t;         ///< typedef to this type
      typedef linear_solver_base<strategy_t>      base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<this_t, true>             sp_this_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<linear_solver_prop, true> sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>           sp_matrix_t;    ///< short name to smart pointer on matrix_t

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~gmres_solver2 ();

      //! solve
      /*virtual int solve (matrix_t *matrix, typename strategy_t::item_array_t & rhs,
                         typename strategy_t::item_array_t & sol);*/

      virtual int solve (matrix_t *matrix, rhs_item_array_t & rhs,
                        item_array_t & sol);

      virtual int solve_prec (matrix_t *matrix, item_array_t & rhs,
                        item_array_t & sol);

      template <class rhs_t>
      int
      templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol);

      //! setup
      virtual int setup (matrix_t *matrix);

    private:
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
      int m;              //!< number of directions to store
      std::vector <item_array_t> vec_p;
      item_array_t vec_w;
      item_array_t vec_r;

    public:
      BLUE_SKY_TYPE_DECL (gmres_solver2<strategy_t>);
    };


  /**
  * @brief BiCGStab linear solver
  */
  template <class strategy_t>
  class BS_API_PLUGIN bicgstab_solver : public linear_solver_base<strategy_t>
    {

      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      typedef typename strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
      typedef typename strategy_t::item_array_t   item_array_t;   ///< short name to array type
      typedef typename strategy_t::item_t         item_t;         ///< short name to array item type
      typedef typename strategy_t::index_t        index_t;        ///< short name to matrix's index type

      typedef typename strategy_t::rhs_item_t	  rhs_item_t;     ///< short name for type of rhs
      typedef typename strategy_t::rhs_item_array_t	rhs_item_array_t;   ///< short name for rhs array type

      typedef typename strategy_t::index_array_t  index_array_t;

      typedef linear_solver_base<strategy_t>      this_t;         ///< typedef to this type
      typedef linear_solver_base<strategy_t>      base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<this_t, true>             sp_this_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<linear_solver_prop, true> sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>           sp_matrix_t;    ///< short name to smart pointer on matrix_t

      typedef bcsr_matrix<rhs_item_array_t, index_array_t>    bcsr_matrix_t;        ///< short name for used matrix
      typedef smart_ptr<bcsr_matrix_t, true>                  sp_bcsr_matrix_t;     ///< short name for smart_pointer on used matrix


      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~bicgstab_solver ();

      //! solve
      //virtual int solve (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);
      virtual int solve (matrix_t *matrix, rhs_item_array_t & rhs,
                         item_array_t & sol);

      virtual int solve_prec (matrix_t *matrix, item_array_t & rhs,
                         item_array_t & sol);

      template <class rhs_t>
      int
      templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol);


      //! setup
      virtual int setup (matrix_t *matrix);

    private:
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:

    public:
      BLUE_SKY_TYPE_DECL (bicgstab_solver);
    };

  //! register linear_solver_prop into kernel
  bool
  linear_solver_prop_register_type (const blue_sky::plugin_descriptor &pd);

  //! register linear_solvers into kernel
  bool
  linear_solvers_register_type (const blue_sky::plugin_descriptor &pd);

}	// namespace blue_sky
#endif //__LIN_SOLVER__H
