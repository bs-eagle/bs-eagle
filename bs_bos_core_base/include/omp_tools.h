#ifndef __OMP_TOOLS_H_
#define __OMP_TOOLS_H_
/*!
 * \file omp_tools.h
 * \brief omp tools and definitions
 * \author Khait Mark
 * \date 2006-07-21
 */

//#define OMP_DEBUG

#ifdef UNIX
#ifndef _OPENMP
#undef OMP_DEBUG
#endif
#endif

#ifdef OMP_DEBUG
#include <omp.h>


//! start OMP time measure
#define OMP_TIME_MEASURE_START(t)                                      \
    t -= omp_get_wtime ()

//! end OMP time measure
#define OMP_TIME_MEASURE_END(t)                                        \
    t += omp_get_wtime ()

#else //OMP_DEBUG
//! start OMP time measure
#define OMP_TIME_MEASURE_START(t)                                      \
    t = 0;

//! end OMP time measure
#define OMP_TIME_MEASURE_END(t)                                        \
    t = 0;
#endif //OMP_DEBUG

#ifdef _OPENMP
#define CSR_MATRIX_VECTOR_PRODUCT_PARALLEL
#if 0
#define SMOOTHER_GS_F_ONLY_INV_PARALLEL
#define SMOOTHER_GS_F_ONLY_PARALLEL
#define SMOOTHER_GS_C_ONLY_INV_PARALLEL
#define SMOOTHER_GS_C_ONLY_PARALLEL
#endif

#define FI_OPERATOR_BLOCK_CONNECTIONS_PARALLEL
#define FI_OPERATOR_CELLS_PARALLEL
#define TRIPLE_MATRIX_PRODUCT_PARALLEL
#define CPR_PRECONDITIONER_SETUP_PARALLEL
#define GMRES_SOLVER2_SOLVE_PARALLEL
#define MV_CALC_LIN_COMB_PARALLEL
#define CPR_PRECONDITIONER_SOLVE_PARALLEL
#define CSR_CALC_LIN_COMB_PARALLEL
#define MV_VECTOR_INNER_PRODUCT_PARALLEL
#define OTHER_NON_IMPORTANT_PARALLEL
#define B_MATRIX_VECTOR_PRODUCT_PARALLEL
#define BUILD_STRENGTH_MATRIX_PARALLEL
#define BUILD_P_PARALLEL
#define BUILD_P_STANDART_PARALLEL
#define COARSE_PMIS_2_PARALLEL
//#define ILU_PREC_PARALLEL
//#define AMG_COARSE_RUGE_PARALLEL
//#define AMG_COARSE_CLJP_PARALLEL
//#define TWO_STAGE_PREC_PARALLEL
#endif //_OPENMP


extern double fi_operator_block_connections_timer;
extern double fi_operator_cells_timer;
extern double cpr_preconditioner_setup_timer;
extern double gmres_solver2_solve_timer;
extern double mv_calc_lin_comb_timer;
extern double cpr_preconditioner_solve_timer;
extern double b_matrix_vector_product_timer;
extern double gmres_setup_timer;
extern double gmres_solve_timer;
extern double ilu_combine_matrix;
extern double ilu_setup_timer;
extern double ilu_solve_timer;

#endif //__OMP_TOOLS_H_
