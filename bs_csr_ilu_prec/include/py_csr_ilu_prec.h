/**
 * \file py_csr_ilu_prec.h
 * \brief Python wrapper for csr_ilu_prec
 * \author Miryanov Sergey
 * \date 18.04.2008
 */

#ifndef PY_CSR_ILU_PREC_H_
#define PY_CSR_ILU_PREC_H_

#ifdef BSPY_EXPORTING_PLUGIN

#include "csr_ilu_prec.h"

namespace blue_sky
  {
  namespace python
    {

    //! export wrappers to python
    void
    py_export_csr_ilu_prec ();

  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

#endif // #ifndef PY_CSR_ILU_PREC_H_
