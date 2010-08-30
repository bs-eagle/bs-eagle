/*!
 * \file bcsr_ilu_prec.cpp
 * \brief implementation of incomplit LU precondition for CSR matrix
 * \author Borschuk Oleg
 * \date 2006-11-07
 */
#include "bs_kernel.h"
#include "bs_assert.h"

//#include "bs_csr_ilu_prec_stdafx.h"
#include "bcsr_ilu_prec.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "strategies.h"
#include "matrix_iface.h"
#include "matrix_macroses.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  /*!
   * \brief constructor
   */
  template <class strategy_t>
  bcsr_ilu_prec<strategy_t>::bcsr_ilu_prec (bs_type_ctor_param /*param*/)
//      : sp_ilu(BS_KERNEL.create_object(bcsr_matrix_t::bs_type()))
      //: lsolver_iface <strategy_t> ()
  {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      lu_matrix = 0; 
      init_prop ();
  }

  template <class strat_t>
  bcsr_ilu_prec<strat_t>::bcsr_ilu_prec(const bcsr_ilu_prec& solver)
      : bs_refcounter (solver) //, lsolver_iface <strat_t> ()
  {
    if (&solver != this)
      *this = solver;
  }

  /*!
   * \brief destructor
   */
  template <class strategy_t>
  bcsr_ilu_prec<strategy_t>::~bcsr_ilu_prec ()
  {
  }

  //! set solver's properties
  template <class strat_t>
  void bcsr_ilu_prec<strat_t>::set_prop(sp_prop_t prop_)
  {
    prop = prop_;

    init_prop ();
  }

  template <class strat_t> void
  bcsr_ilu_prec<strat_t>::init_prop ()
    {
      use_internal_matrix_idx = prop->get_index_b (std::string ("use_internal_matrix"));
      if (use_internal_matrix_idx < 0)
        use_internal_matrix_idx = prop->add_property_b (true, std::string ("use_internal_matrix"), 
                                                        std::string ("If True preconditioner will copy matrix given in setup method"));

      if (use_internal_matrix_idx < 0)
        {
          bs_throw_exception ("Can not regidter some properties");
        }
    }

  template <class strategy_t>
  int bcsr_ilu_prec<strategy_t>::solve(sp_matrix_t matrix, sp_fp_array_t sp_rhs, sp_fp_array_t sp_sol)
  {
    BS_ASSERT (matrix);
    BS_ASSERT (sp_rhs->size ());
    BS_ASSERT (sp_sol->size ());
    BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());

    i_type_t b_sqr; 
    i_type_t n;     
    i_type_t nb;    
    fp_type_t *rhs = &(*sp_rhs)[0];
    fp_type_t *sol = &(*sp_sol)[0];

    sp_bcsr_matrix_t ilu;

    if (dynamic_cast<bcsr_matrix_iface_t*> (matrix.lock ()))
      {
        ilu = matrix;
        BS_ASSERT (ilu);
      }

    //const item_array_t &values    = ilu.get_values   ();
    bool ff = prop->get_b (use_internal_matrix_idx);

    sp_fp_storage_array_t sp_values     = ff ? lu_matrix->get_values () : ilu->get_values   ();
    sp_i_array_t sp_rows                = ff ? lu_matrix->get_rows_ptr () : ilu->get_rows_ptr ();           
    sp_i_array_t sp_cols                = ff ? lu_matrix->get_cols_ind () : ilu->get_cols_ind ();           
    fp_storage_type_t *values           = &(*sp_values)[0];
    i_type_t *rows                      = &(*sp_rows)[0];
    i_type_t *cols                      = &(*sp_cols)[0];

    if (ff)
      {
        n         = lu_matrix->get_n_rows ();
        nb        = lu_matrix->get_n_block_size ();
      }
    else
      {

        n         = ilu->get_n_rows ();
        nb        = ilu->get_n_block_size ();
      }
    b_sqr = nb * nb;

    i_type_t j, k, l;
    i_type_t j1, j2;
    //const fp_type *D_block, *M_block;
    
    const fp_storage_type_t *D_block, *M_block;

    fp_type_t *v, *r;

    memcpy (sol, rhs, sp_rhs->size () * sizeof (rhs[0]));

#ifdef _DEBUG
    fp_type_t *sol_ = &sol[0];
#endif
    // solve Ly = b
    r = &sol[0];
    for (k = 0; k < n; ++k, r += nb)
      {
        // pointer to matrix row
        j1          = rows[k];
        j2          = j1;//diag_ind[k];
        i_type_t j3  = rows[k + 1];
        for (j = j1 + 1; j < j3; ++j)
          {
            l = cols[j];
            if (l < k)
              {
                v = &sol[l * nb];
                M_block = &values[j * b_sqr];
                LU_FIND_ROOT_UPDATE (nb, M_block, v, r);
              }
            else
              break;
          }
        BS_ASSERT (j1 == j2);
        BS_ASSERT (cols[j2] == k) (cols[j2]) (k);
        D_block = &values[j2 * b_sqr];
        uLU_FIND_ROOT_L (nb, D_block, r);
      }
    // solve Ux = y
    r = &sol [(n - 1) * nb];
    for (k = n - 1; k >= 0; --k, r -= nb)
      {
        // find element in k column
        i_type_t j0  = rows[k];
        j1          = j0;//diag_ind[k];
        j2          = rows[k + 1];
        for (j = j0; j < j2; ++j)
          {
            l = cols[j];
            if (l > k)
              {
                v = &sol[l * nb];
                M_block = &values[j * b_sqr];
                LU_FIND_ROOT_UPDATE (nb, M_block, v, r);
              }
          }
        BS_ASSERT (j0 == j1);
        BS_ASSERT (cols[j1] == k) (cols[j1]) (k);
        D_block = &values[j1 * b_sqr];
        uLU_FIND_ROOT_U (nb, D_block, r);
      }

    return 0;
  }

  template <class strategy_t>
  int bcsr_ilu_prec<strategy_t>::solve_prec(sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol)
  {
     return solve (matrix, rhs, sol);
  }

  /**
   * \brief setup preconditioner (merge matrices if needed)
   *
   * \param matrix Various matrix
   * \return 0 if success
   */
  template <class strategy_t> int
  bcsr_ilu_prec<strategy_t>::setup (sp_matrix_t matrix)
  {
    BS_ASSERT (matrix);


    i_type_t n;     
    i_type_t nb;    
    i_type_t b_sqr; 
    bool ff = prop->get_b (use_internal_matrix_idx);


    sp_bcsr_matrix_t ilu;

    if (dynamic_cast<bcsr_matrix_iface_t *> (matrix.lock ()))
      {
        ilu = matrix;
        BS_ASSERT (ilu);
      }

    if (ff)
      {
        if (!lu_matrix)
          {
            lu_matrix = BS_KERNEL.create_object (matrix->bs_resolve_type ());
            if (!lu_matrix)
              {
                bs_throw_exception ("Can not create matrix");
              }
          }
        if (lu_matrix->copy (ilu))
          {
            bs_throw_exception ("Can not make matrix copy");
          }

        n              = lu_matrix->get_n_rows ();
        nb             = lu_matrix->get_n_block_size ();
      }
    else
      {
        n              = ilu->get_n_rows ();
        nb             = ilu->get_n_block_size ();
      }

    sp_i_array_t sp_ilu_rows                    = ff ? lu_matrix->get_rows_ptr () : ilu->get_rows_ptr ();       
    sp_i_array_t sp_ilu_cols                    = ff ? lu_matrix->get_cols_ind () : ilu->get_cols_ind ();       
    sp_fp_storage_array_t  sp_ilu_values        = ff ? lu_matrix->get_values ()   : ilu->get_values (); 

    i_type_t *ilu_rows                          = &(*sp_ilu_rows)[0];
    i_type_t *ilu_cols                          = &(*sp_ilu_cols)[0];
    fp_storage_type_t *ilu_values               = &(*sp_ilu_values)[0];

    b_sqr          = nb * nb;

    // common
    //int r_code = 0;
    i_type_t i, j, j1, j2, j3, cl, k;

    fp_storage_type_t *block;



    //ilu.ascii_write_in_csr_format ("ilu_set_internal");

    //  ilu.write_matrix_to_file ("matrix.out");
    // ----------------------
    // STAGE 3: build ilu factorization
    i_type_t i_str;
    //fp_type *d_block, *dd_block;
    fp_storage_type_t *d_block, *dd_block;
    fp_type_t d;
    i_type_t jj, jj1, jj2;

    // loop through all rows
    for (i = 0; i < n; ++i)
      {
        // update i-th row by 0..i-1 row
        j1 = ilu_rows[i];
        j2 = j1;//ilu_diag_ind[i];
        j3 = ilu_rows[i + 1];

        BS_ASSERT (j1 == j2) (j1) (j2);

        // for all nonzero elements in i-th string (in L part)
        // update it by diagonal element and update all other nonzero elements
        // in i-th string
        for (j = j1 + 1; j < j3; ++j)
          {
            // find corresponding diagonal element
            i_str = ilu_cols[j];
            if (i_str < i)
              {
                // update by diagonal element
                block   = &ilu_values[j * b_sqr];
                d_block = &ilu_values[ilu_rows[i_str]/*ilu_diag_ind[i_str]*/ * b_sqr];
                //BS_ASSERT (ilu_diag_ind[i_str] == ilu_rows[i_str]) (ilu_diag_ind[i_str]) (ilu_rows[i_str]);
                uLU_SEEK_L (nb, d_block, block, d);

                jj1 = ilu_rows[i_str];
                jj2 = ilu_rows[i_str + 1];
                k = j + 1;

                for (jj = jj1; jj < jj2; ++jj)
                  {
                    cl = ilu_cols[jj];
                    if (cl <= i_str)
                      continue;

                    if (ilu_cols[j1] == cl)
                      {
                        // upgrade by corresponding element
                        d_block   = &ilu_values[jj * b_sqr];
                        dd_block  = &ilu_values[j1 * b_sqr];
                        LU_UPGRADE(nb, block, d_block, dd_block);
                      }
                    for (; k < j3 && ilu_cols[k] < cl; ++k)
                      ;

                    if (k < j3 && ilu_cols[k] == cl)
                      {
                        // upgrade by corresponding element
                        d_block   = &ilu_values[jj * b_sqr];
                        dd_block  = &ilu_values[k * b_sqr];
                        LU_UPGRADE(nb, block, d_block, dd_block);
                      }
                    else if (k >= j3)
                      break;
                  }
              }
            else 
              break;
          }

        // factorize i-th string
        block = &ilu_values[j2 * b_sqr];
        uLU (nb, block, d);

        for (j = j1; j < j3; ++j)
          {
            i_type_t cl = ilu_cols[j];
            if (cl > i)
              {
                d_block = &ilu_values[j * b_sqr];
                uLU_SEEK_U (nb, block, d_block);
              }
          }
      }
    // ----------------------
    return 0;
  }


  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(bcsr_ilu_prec, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(bcsr_ilu_prec, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_ilu_prec<base_strategy_fif>) , 1, (lsolver_iface<base_strategy_fif>), "bcsr_ilu_prec_fif", "BCSR ILU Preconditioner", "BCSR ILU Preconditioner", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_ilu_prec<base_strategy_did>) , 1, (lsolver_iface<base_strategy_did>), "bcsr_ilu_prec_did", "BCSR ILU Preconditioner", "BCSR ILU Preconditioner", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_ilu_prec<base_strategy_dif>) , 1, (lsolver_iface<base_strategy_dif>), "bcsr_ilu_prec_dif", "BCSR ILU Preconditioner", "BCSR ILU Preconditioner", false);

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_ilu_prec<base_strategy_flf>) , 1, (lsolver_iface<base_strategy_flf>), "bcsr_ilu_prec_flf", "BCSR ILU Preconditioner", "BCSR ILU Preconditioner", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_ilu_prec<base_strategy_dld>) , 1, (lsolver_iface<base_strategy_dld>), "bcsr_ilu_prec_dld", "BCSR ILU Preconditioner", "BCSR ILU Preconditioner", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bcsr_ilu_prec<base_strategy_dlf>) , 1, (lsolver_iface<base_strategy_dlf>), "bcsr_ilu_prec_dlf", "BCSR ILU Preconditioner", "BCSR ILU Preconditioner", false);

#if 0
  //////////////////////////////////////////////////////////////////////////
  //! register types in kernel
  bool bcsr_ilu_prec_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, bcsr_ilu_prec<base_strategy_fi>::bs_type ());
    res &= BS_KERNEL.register_type (pd, bcsr_ilu_prec<base_strategy_di>::bs_type ());
    res &= BS_KERNEL.register_type (pd, bcsr_ilu_prec<base_strategy_mixi>::bs_type ());



    return res;
  }
#endif //0

} // namespace blue_sky
