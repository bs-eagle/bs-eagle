/**
 * \file   setup_preconditioner.h
 * \brief  Base interface for setup various matrices for using in preconditioner
 * \author Miryanov Sergey
 * \date 2008-04-14
 */

#ifndef BS_SETUP_PRECONDITIONER_H_
#define BS_SETUP_PRECONDITIONER_H_

#include "matrix_base.h"

namespace blue_sky
  {


  /**
   * \brief setup_preconditioner
   */
  template <class matrix_t>
  class setup_preconditioner
    {
    private:

      typedef smart_ptr<matrix_t, true> sp_matrix_t;  ///< short name for smart pointer on matrix

    public:

      /**
       * \brief copy current matrix and prepare it for using in any preconditioners
       *
       * \return The BCSR matrix or 0
       */
      virtual void prepare_matrix () = 0;

      virtual sp_matrix_t get_prepared_matrix () const = 0;
    };

}


#endif // #ifndef BS_SETUP_PRECONDITIONER_H_
