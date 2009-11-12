/**
 *       \file  fi_operator_norm_calc.h
 *      \brief  Calculates norm for fi_operator
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.01.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_FI_OPERATOR_NORM_CALC_H_
#define BS_FI_OPERATOR_NORM_CALC_H_

namespace blue_sky {

#define MAX_AND_INDEX(NEW,OLD,I,I_OLD) if (fabs (NEW) > fabs (OLD)) {(OLD) = (NEW); (I_OLD) = (I);}
#ifdef _MPI
struct di
  {
    double d;
    int i;
  };
#define MPI_MAX_AND_INDEX(NEW,I,DI_OLD,DI_NEW)  DI_OLD.d = fabs (NEW);  DI_OLD.i = (NEW > 0) ? I : -I;\
  MPI_Allreduce (&DI_OLD, &DI_NEW, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);  \
  NEW = (DI_NEW.i > 0) ? DI_NEW.d : - DI_NEW.d; I = abs(DI_NEW.i);
#endif //_MPI

  /**
   * \brief Calculates norms
   */
  template <typename strategy_t, bool is_w, bool is_g, bool is_o>
  BS_FORCE_INLINE void
  fi_operator_impl <strategy_t, is_w, is_g, is_o>::norm_calc ()
  {
    double pv = 0;

#ifdef _MPI
    int n_left = 0;///mpi_decomp->get_recv_l ();
    int n_right = n_cells_;///mpi_decomp->get_n_local_own () + n_left;
    int n_offset = 0;///mpi_decomp->get_elem_own_displ () + n_left;
    double mpi_l2_norm[6], mpi_mat_balance[3];
    struct di max_ind_old, max_ind_new;
#else //_MPI
    int n_left = 0;
    int n_right = n_cells_;
#endif //_MPI

    norm_.clear ();
    for (index_t i = n_left; i < n_right; ++i)
      {
        // calculate total poro volume
        pv += volume_[i] * data_[i].porosity;
        update_norm_by_cell (i, norm_);
      }

#ifdef _MPI
    BS_ASSERT (false && "MPI: NOT IMPL YET");
    // get global pv
    double mpi_pv;
    MPI_Allreduce (&pv, &mpi_pv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    pv = mpi_pv;

    ///TODO: MPI_TYPE check!

    // collecting l2 norms from all procs
    MPI_Allreduce (&norm.val[norms::L2_CPV_WATER], &mpi_l2_norm[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::L2_CPV_GAS], &mpi_l2_norm[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::L2_CPV_OIL], &mpi_l2_norm[2], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::L2_ACPV_WATER], &mpi_l2_norm[3], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::L2_ACPV_GAS], &mpi_l2_norm[4], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::L2_ACPV_OIL], &mpi_l2_norm[5], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::MB_ERR_WATER], &mpi_mat_balance[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::MB_ERR_GAS], &mpi_mat_balance[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&norm.val[norms::MB_ERR_OIL], &mpi_mat_balance[2], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    norm.val[norms::L2_CPV_WATER] = mpi_l2_norm[0];
    norm.val[norms::L2_CPV_GAS] = mpi_l2_norm[1];
    norm.val[norms::L2_CPV_OIL] = mpi_l2_norm[2];
    norm.val[norms::L2_ACPV_WATER] = mpi_l2_norm[3];
    norm.val[norms::L2_ACPV_GAS] = mpi_l2_norm[4];
    norm.val[norms::L2_ACPV_OIL] = mpi_l2_norm[5];
    norm.val[norms::MB_ERR_WATER] = mpi_mat_balance[0];
    norm.val[norms::MB_ERR_GAS] = mpi_mat_balance[1];
    norm.val[norms::MB_ERR_OIL] = mpi_mat_balance[2];

    MPI_MAX_AND_INDEX (norm.val[norms::C_CPV_WATER], norm.idx[norms::C_CPV_WATER], max_ind_old, max_ind_new);
    MPI_MAX_AND_INDEX (norm.val[norms::C_CPV_GAS], norm.idx[norms::C_CPV_GAS], max_ind_old, max_ind_new);
    MPI_MAX_AND_INDEX (norm.val[norms::C_CPV_OIL], norm.idx[norms::C_CPV_OIL], max_ind_old, max_ind_new);
    MPI_MAX_AND_INDEX (norm.val[norms::C_ACPV_WATER], norm.idx[norms::C_ACPV_WATER], max_ind_old, max_ind_new);
    MPI_MAX_AND_INDEX (norm.val[norms::C_ACPV_GAS], norm.idx[norms::C_ACPV_GAS], max_ind_old, max_ind_new);
    MPI_MAX_AND_INDEX (norm.val[norms::C_ACPV_OIL], norm.idx[norms::C_ACPV_OIL], max_ind_old, max_ind_new);

    ///n = mpi_decomp->get_n_elements ();
#endif //_MPI

    double n = 1.0 / double (n_cells_);
    norm_.val[norms::L2_CPV_WATER]   = sqrt (norm_.val[norms::L2_CPV_WATER])  * n;
    norm_.val[norms::L2_CPV_GAS]     = sqrt (norm_.val[norms::L2_CPV_GAS])    * n;
    norm_.val[norms::L2_CPV_OIL]     = sqrt (norm_.val[norms::L2_CPV_OIL])    * n;
    norm_.val[norms::L2_ACPV_WATER]  = sqrt (norm_.val[norms::L2_ACPV_WATER]) * n;
    norm_.val[norms::L2_ACPV_GAS]    = sqrt (norm_.val[norms::L2_ACPV_GAS])   * n;
    norm_.val[norms::L2_ACPV_OIL]    = sqrt (norm_.val[norms::L2_ACPV_OIL])   * n;

    // C_CPV and
    if (is_w)
      {
        MAX_AND_INDEX (norm_.val[norms::C_CPV_WATER],   norm_.val[norms::C_CPV],
                       norm_.idx[norms::C_CPV_WATER],         norm_.idx[norms::C_CPV]);
        norm_.val[norms::C_ACPV_WATER] /= calc_model_->ave_volume;

        MAX_AND_INDEX (norm_.val[norms::C_ACPV_WATER],  norm_.val[norms::C_ACPV],
                       norm_.idx[norms::C_ACPV_WATER],        norm_.idx[norms::C_ACPV]);

        MAX_AND_INDEX (norm_.val[norms::L2_CPV_WATER],  norm_.val[norms::L2_CPV],
                       norm_.idx[norms::L2_CPV_WATER],        norm_.idx[norms::L2_CPV]);

        MAX_AND_INDEX (norm_.val[norms::L2_ACPV_WATER], norm_.val[norms::L2_ACPV],
                       norm_.idx[norms::L2_ACPV_WATER],       norm_.idx[norms::L2_ACPV]);

        norm_.val[norms::MB_ERR_WATER] /= pv * calc_model_->invers_fvf_average[d_w];
        MAX_AND_INDEX (norm_.val[norms::MB_ERR_WATER],  norm_.val[norms::MB_ERR],
                       norm_.idx[norms::MB_ERR_WATER],        norm_.idx[norms::MB_ERR]);
      }
    if (is_o)
      {
        MAX_AND_INDEX (norm_.val[norms::C_CPV_OIL],     norm_.val[norms::C_CPV],
                       norm_.idx[norms::C_CPV_OIL],           norm_.idx[norms::C_CPV]);
        norm_.val[norms::C_ACPV_OIL] /= calc_model_->ave_volume;

        MAX_AND_INDEX (norm_.val[norms::C_ACPV_OIL],    norm_.val[norms::C_ACPV],
                       norm_.idx[norms::C_ACPV_OIL],          norm_.idx[norms::C_ACPV]);

        MAX_AND_INDEX (norm_.val[norms::L2_CPV_OIL],    norm_.val[norms::L2_CPV],
                       norm_.idx[norms::L2_CPV_OIL],          norm_.idx[norms::L2_CPV]);

        MAX_AND_INDEX (norm_.val[norms::L2_ACPV_OIL],   norm_.val[norms::L2_ACPV],
                       norm_.idx[norms::L2_ACPV_OIL],         norm_.idx[norms::L2_ACPV]);

        norm_.val[norms::MB_ERR_OIL] /= pv * calc_model_->invers_fvf_average[d_o];
        MAX_AND_INDEX (norm_.val[norms::MB_ERR_OIL],    norm_.val[norms::MB_ERR],
                       norm_.idx[norms::MB_ERR_OIL],          norm_.idx[norms::MB_ERR]);

      }
    if (is_g)
      {
        MAX_AND_INDEX (norm_.val[norms::C_CPV_GAS],     norm_.val[norms::C_CPV],
                       norm_.idx[norms::C_CPV_GAS],           norm_.idx[norms::C_CPV]);
        norm_.val[norms::C_ACPV_GAS] /= calc_model_->ave_volume;

        MAX_AND_INDEX (norm_.val[norms::C_ACPV_GAS],    norm_.val[norms::C_ACPV],
                       norm_.idx[norms::C_ACPV_GAS],          norm_.idx[norms::C_ACPV]);

        MAX_AND_INDEX (norm_.val[norms::L2_CPV_GAS],    norm_.val[norms::L2_CPV],
                       norm_.idx[norms::L2_CPV_GAS],          norm_.idx[norms::L2_CPV]);

        MAX_AND_INDEX (norm_.val[norms::L2_ACPV_GAS],   norm_.val[norms::L2_ACPV],
                       norm_.idx[norms::L2_ACPV_GAS],         norm_.idx[norms::L2_ACPV]);

        norm_.val[norms::MB_ERR_GAS] /= pv * calc_model_->invers_fvf_average[d_g];
        MAX_AND_INDEX (norm_.val[norms::MB_ERR_GAS],    norm_.val[norms::MB_ERR],
                       norm_.idx[norms::MB_ERR_GAS],          norm_.idx[norms::MB_ERR]);
      }


#ifdef _MPI
    if ((norm_.idx[norms::C_CPV] - n_offset >= n_left) && (norm_.idx[norms::C_CPV] - n_offset < n_right))
      {
        calc_model_->max_norm_counter[norm_.idx[norms::C_CPV] - n_offset]++;
      }
#else //_MPI
    calc_model_->max_norm_counter[norm_.idx[norms::C_CPV]]++;
#endif //_MPI
  }

  /**
   * \brief  Calculates norm in one cell and update norm storage
   * \param  i Index of cell
   * \param  ns Norms storage
   * */
  template <class strategy_t, bool is_w, bool is_g, bool is_o>
  BS_FORCE_INLINE void 
  fi_operator_impl <strategy_t, is_w, is_g, is_o>::update_norm_by_cell (index_t i, norms_storage_t &ns)
  {
    double norm_mult, r, d;

    int equ_w, equ_g, equ_o;
    if (n_phases == 3)
      {
        equ_w = i * n_phases + p3_wat;
        equ_g = i * n_phases + p3_gas;
        equ_o = i * n_phases + p3_oil;
      }
    else if (n_phases == 2 && is_w)
      {
        equ_w = i * n_phases + p2ow_wat;
        equ_o = i * n_phases + p2ow_oil;
      }
    else if (n_phases == 2 && is_g)
      {
        equ_g = i * n_phases + p2og_gas;
        equ_o = i * n_phases + p2og_oil;
      }
    else
      {
        equ_g = equ_w = equ_o = i * n_phases + 0;
      }
#ifdef _MPI
    int n_left = 0;///mpi_decomp->get_recv_l ();
    int n_offset = 0;///mpi_decomp->get_elem_own_displ () + n_left;
    ///rhs -= n_phases * n_left;
    ///flux_rhs -= n_phases * n_left;
#endif //_MPI

    // calculate multiplier
    norm_mult = 1.0 / (volume_[i] * data_[i].porosity);

#ifdef _MPI
    // remember absolute (global) cell number for norms
    i += n_offset;
#endif //_MPI

    // water
    if (is_w)
      {
        r = rhs_[equ_w] + flux_rhs_[equ_w];

        // mat balanse
        ns.val[norms::MB_ERR_WATER] += r / calc_model_->invers_fvf_average[d_w];

        // C_CPV_WATER
        d = r * norm_mult / calc_model_->invers_fvf_average[d_w];
        MAX_AND_INDEX (d, ns.val[norms::C_CPV_WATER], i, ns.idx[norms::C_CPV_WATER]);

        // L2_CPV_WATER
        ns.val[norms::L2_CPV_WATER] += d * d;

        // C_ACPV_WATER
        d = r / data_[i].invers_fvf[d_w];
        MAX_AND_INDEX (d, ns.val[norms::C_ACPV_WATER], i, ns.idx[norms::C_ACPV_WATER]);

        // L2_ACPV_WATER
        ns.val[norms::L2_ACPV_WATER] += d * d;
      }
    // gas
    if (is_g)
      {
        r = rhs_[equ_g] + flux_rhs_[equ_g];

        // mat balanse
        ns.val[norms::MB_ERR_GAS] += r / calc_model_->invers_fvf_average[d_g];

        // C_CPV_GAS
        d = r * norm_mult / calc_model_->invers_fvf_average[d_g];
        MAX_AND_INDEX (d, ns.val[norms::C_CPV_GAS], i, ns.idx[norms::C_CPV_GAS]);

        // L2_CPV_GAS
        ns.val[norms::L2_CPV_GAS] += d * d;

        // C_ACPV_GAS
        d = r / data_[i].invers_fvf[d_g];
        MAX_AND_INDEX (d, ns.val[norms::C_ACPV_GAS], i, ns.idx[norms::C_ACPV_GAS]);

        // L2_ACPV_GAS
        ns.val[norms::L2_ACPV_GAS] += d * d;
      }
    // oil
    if (is_o)
      {
        r = rhs_[equ_o] + flux_rhs_[equ_o];

        // mat balanse
        ns.val[norms::MB_ERR_OIL] += r / calc_model_->invers_fvf_average[d_o];

        // C_CPV_OIL
        d = r * norm_mult / calc_model_->invers_fvf_average[d_o];
        MAX_AND_INDEX (d, ns.val[norms::C_CPV_OIL], i, ns.idx[norms::C_CPV_OIL]);

        // L2_CPV_OIL
        ns.val[norms::L2_CPV_OIL] += d * d;

        // C_ACPV_OIL
        d = r / data_[i].invers_fvf[d_o];
        MAX_AND_INDEX (d, ns.val[norms::C_ACPV_OIL], i, ns.idx[norms::C_ACPV_OIL]);

        // L2_ACPV_OIL
        ns.val[norms::L2_ACPV_OIL] += d * d;
      }

  rhs_item_array_t &s_rhs = jmatrix_->get_sec_rhs ();
  for (index_t j = 0, j_cnt = n_sec_vars; j < j_cnt; ++j)
    {
      MAX_AND_INDEX (s_rhs[i * n_sec_vars + j], ns.val[norms::S_RHS], i, ns.idx[norms::S_RHS]);
    }
  }

} // namespace blue_sky

#endif // #ifndef BS_FI_OPERATOR_NORM_CALC_H_
