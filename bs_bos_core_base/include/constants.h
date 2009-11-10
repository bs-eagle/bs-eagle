/**
 * \file constants.h
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */
#ifndef BS_SOLVER_CONSTANTS_H_
#define BS_SOLVER_CONSTANTS_H_

namespace blue_sky
  {

#define EPS_DIFF  1.0e-7
#define EPS_DIV   1e-15


  /*!
    \brief variables in grid block for gas-phase model
  */

  enum well_model_type
  {
    BLACK_OIL = 0,
    COMPOSIT
  };

  enum main_var_type
  {
    FI_NULL = -1,
    FI_SG_VAR = 0,
    FI_RO_VAR,
    FI_MOMG_VAR,
  };

  //! well model type
  enum well_model_var_type
  {
    WELL_MODEL_1VAR = 0,
    WELL_MODEL_3VAR
  };

  enum RPO_MODEL_ENUM
  {
    RPO_DEFAULT_MODEL = 0,
    STONE1_MODEL,
    STONE2_MODEL
  };

  enum PHASE_ENUM
  {
    PHASE_NULL		= 0,
    PHASE_WATER   = 1,
    PHASE_GAS     = 2,
    PHASE_OIL     = 4,

    PHASE_TOTAL,
  };

  enum FI_PHASE_ENUM
  {
    FI_PHASE_NULL		= 0,
    FI_PHASE_WATER  = 0,
    FI_PHASE_GAS,
    FI_PHASE_OIL,
    FI_PHASE_TOT,
  };

  //! types of gas phase rates
  enum
  {
    GAS_RATE_FREE,
    GAS_RATE_SOL,
    GAS_RATE_TOTAL
  };

  enum
  {
    FI_WELL_FW_VAR = 0,
    FI_WELL_FG_VAR,
    FI_WELL_XREF_VAR
  };

  enum FI_LIN_SOLVER_ENUM
  {
    FI_LIN_SOLVER_BICGSTAB = 0,
    FI_LIN_SOLVER_GMRES
  };

  enum FI_LIN_PREC_ENUM
  {
    FI_LIN_PREC_CPR = 0,
    FI_LIN_PREC_ILU,
    FI_LIN_PREC_CPR_SOR,
    FI_LIN_PREC_AMG
  };

  // check if phase is present
#ifdef COMPOSITIONAL
#define FI_CHK_WATER(PH)        ((PH) & (1 << COMP_PHASE_WATER))
#define FI_CHK_GAS(PH)          ((PH) & (1 << COMP_PHASE_GAS))
#define FI_CHK_OIL(PH)          ((PH) & (1 << COMP_PHASE_OIL))
#else
#define FI_CHK_WATER(PH)        ((PH) & (1 << FI_PHASE_WATER))
#define FI_CHK_GAS(PH)          ((PH) & (1 << FI_PHASE_GAS))
#define FI_CHK_OIL(PH)          ((PH) & (1 << FI_PHASE_OIL))
#endif

  // get index I of phase with displacement D in array, N -- total number of phases
#define FI_PH_IND(I,D,N)        ((I) * (N) + (D))

  // check for oil - gas gas system is present
#define FI_CHK_OIL_GAS(PH)      (FI_CHK_GAS((PH)) && FI_CHK_OIL((PH)))
  // check for oil - water system is present
#define FI_CHK_OIL_WATER(PH)    (FI_CHK_OIL((PH)) && FI_CHK_WATER((PH)))
  // check for gas - water system is present
#define FI_CHK_GAS_WATER(PH)    (FI_CHK_GAS((PH)) && FI_CHK_WATER((PH)))

  //---------------------
  //  3 phase (gas, oil, water) system
  //---------------------
  //! displacements of d()/d(Sg) in Jacobian block row
  const int p3_sg = 0;
  //! displacements of d()/d(So) in Jacobian block row
  const int p3_so = 1;
  //! displacements of d()/d(Po) in Jacobian block row
  const int p3_po = 2;

  //! displacements of d(Rg)/d() in Jacobian block column
  const int p3_gas = 0;
  //! displacements of d(Ro)/d() in Jacobian block column
  const int p3_oil = 1;
  //! displacements of d(Rw)/d() in Jacobian block column
  const int p3_wat = 2;

  enum p3_deriv_index_type
  {
    p3_gas_sg = 0,
    p3_gas_so,
    p3_gas_po,
    p3_oil_sg,
    p3_oil_so,
    p3_oil_po,
    p3_wat_sg,
    p3_wat_so,
    p3_wat_po
  };
  //==============================

  //---------------------
  //  2 phase (gas, oil) system
  //---------------------
  //! displacements of d()/d(Sg) in Jacobian block row
  const int p2og_sg = 0;
  //! displacements of d()/d(Po) in Jacobian block row
  const int p2og_po = 1;

  //! displacements of d(Rg)/d() in Jacobian block column
  const int p2og_gas = 0;
  //! displacements of d(Ro)/d() in Jacobian block column
  const int p2og_oil = 1;

  enum
  {
    p2og_gas_sg = 0,
    p2og_gas_po,
    p2og_oil_sg,
    p2og_oil_po
  };
  //==============================

  //---------------------
  //  2 phase (oil, water) system
  //---------------------
  //! displacements of d()/d(So) in Jacobian block row
  const int p2ow_so = 0;
  //! displacements of d()/d(Po) in Jacobian block row
  const int p2ow_po = 1;

  //! displacements of d(Ro)/d() in Jacobian block column
  const int p2ow_oil = 0;
  //! displacements of d(Rw)/d() in Jacobian block column
  const int p2ow_wat = 1;

  //! values returned by fi_operator method
  enum fi_operator_return_type
  {
    FI_OPERATOR_RETURN_FAIL = 0,
    FI_OPERATOR_RETURN_OK = 1,
    FI_OPERATOR_RETURN_RESTART,
    FI_OPERATOR_APPROX_DT
  };

  enum
  {
    p2ow_oil_so = 0,
    p2ow_oil_po,
    p2ow_wat_so,
    p2ow_wat_po
  };
  //==============================

  //! check for main variable (Sg or Ro) in grid block
#define FI_GET_VAR(VAR,I)       ((VAR)[(I)])
#define FI_CHK_SG(VAR,I)        ((VAR)[(I)]) == FI_SG_VAR ? 1 : 0
#define FI_CHK_RO(VAR,I)        ((VAR)[(I)]) == FI_RO_VAR ? 1 : 0
#define FI_CHK_MOMG(VAR,I)      ((VAR)[(I)]) == FI_MOMG_VAR ? 1 : 0

//! EQUIL keyword enumerated params
  enum
  {
    EQUIL_DAT_DEPTH = 0,
    EQUIL_DAT_PRESS,
    EQUIL_WOC_DEPTH,
    EQUIL_WOC_PRESS,
    EQUIL_GOC_DEPTH,
    EQUIL_GOC_PRESS,
    EQUIL_RS_TYPE,
    EQUIL_RV_TYPE,
    EQUIL_NUM_SEC,
    EQUIL_COMP_TYPE,
    EQUIL_COMP_ARG,
    EQUIL_TOTAL
  };

} // namespace blue_sky


#endif  // #ifndef BS_SOLVER_CONSTANTS_H_
