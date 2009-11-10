#include "bs_matrix_stdafx.h"

#ifdef _MPI
#include "mpi_vector.h"
#endif

#include "matrix_base.h"
#include "bs_kernel.h"
#include "seq_vector.h"

namespace blue_sky
  {


  //! constructor
  template< class fp_vector_type, class i_vector_type>
  matrix_base< fp_vector_type, i_vector_type >::matrix_base(bs_type_ctor_param)
  {
    n_block_size = 1;
    n_rows = 0;
    n_cols = 0;
  }

  template< class fp_vector_type, class i_vector_type >
  matrix_base< fp_vector_type, i_vector_type >::matrix_base(const matrix_base & src) 
        : bs_refcounter (), objbase ()
  {
    *this = src;
  }

  //Blue sky register stuff
  //create instances

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(matrix_base, (class)(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(matrix_base, (class)(class));

  BLUE_SKY_TYPE_IMPL_T_EXT(2, (matrix_base<seq_vector<float>, seq_vector<int> >) , 1, (objbase), "matrix_base<float, int>", "Base Matrix class", "Realization of Base Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(2, (matrix_base<seq_vector<double>, seq_vector<int> >) , 1, (objbase), "matrix_base<double, int>", "Base Matrix class", "Realization of Base Matricies", false);

#ifdef _MPI
  BLUE_SKY_TYPE_IMPL_T_EXT(2, (matrix_base<mpi_vector<float>, mpi_vector<int> >) , 1, (objbase), "matrix_mpi<float, int>", "MPI Matrix class", "Realization of MPI Matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(2, (matrix_base<mpi_vector<double>, mpi_vector<int> >) , 1, (objbase), "matrix_mpi<double, int>", "MPI Matrix class", "Realization of MPI Matricies", false);
#endif

}//ns bs

