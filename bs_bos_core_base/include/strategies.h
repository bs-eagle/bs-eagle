/**
* \file   call_helper.cpp
* \brief  Strategies for linear solvers
* \author Miryanov Sergey
* \date 2008-04-03
*/

#ifndef BS_LINEAR_SOLVERS_MATRIX_STRATEGIES_H_
#define BS_LINEAR_SOLVERS_MATRIX_STRATEGIES_H_

#ifdef _MPI
#include <mpi.h>
#include "mpi_vector.h"
#endif

#include "shared_vector.h"

namespace blue_sky
{
  template <typename array_item_t, typename array_index_t>
  class BS_API_PLUGIN matrix_base;

  template <typename array_item_t, typename array_index_t>
  class BS_API_PLUGIN bcsr_matrix;

  template <typename array_item_t, typename array_index_t>
  class BS_API_PLUGIN mpi_csr_matrix;

  struct BS_API_PLUGIN seq_barrier
  {
    typedef int mpi_comm_t;

    seq_barrier ()
    {
    }
    void barrier (mpi_comm_t /* not used */)
    {
    }

    static mpi_comm_t
    comm_world_v ()
    {
      return 0;
    }
  };


#ifdef _MPI
  struct BS_API_PLUGIN mpi_barrier
  {
    typedef MPI_Comm mpi_comm_t;

    mpi_barrier ()
    {
    }
    void barrier (mpi_comm_t comm)
    {
      MPI_Barrier (comm);
    }

    static mpi_comm_t
    comm_world_v ()
    {
      return MPI_COMM_WORLD;
    }
  };
#endif

  template <class rhs_fp_type, class fp_type, class i_type>
  struct BS_API_PLUGIN base_strategy
  {
    typedef seq_barrier                                       barrier_t;          ///< short name for barrier

    typedef fp_type                                           item_t;             ///< short name for floating point type (old array_item_t)
    typedef i_type                                            index_t;            ///< short name for index type
    typedef shared_vector <item_t>                            item_array_t;       ///< short name for array type (old array_t)
    typedef shared_vector <index_t>                           index_array_t;      ///< short name for type of array of indexes
    typedef rhs_fp_type		                                    rhs_item_t;         ///< short name for type of rhs
    typedef shared_vector <rhs_item_t>	                      rhs_item_array_t;   ///< short name for rhs array type

    typedef matrix_base<rhs_item_array_t, index_array_t>      matrix_t;           ///< short name for matrix type
    typedef bcsr_matrix<rhs_item_array_t, index_array_t>      csr_matrix_t;       ///< short name for csr_matrix type (it's fail way today)

    typedef shared_vector <fp_type>                           wksp_t;             ///< short name for workspace array

    //typedef std::vector<matrix_t>                             matrix_array_t;     ///< short name for vector of matrices type

    template <typename T>
    struct vec
      {
        typedef shared_vector <T>                             type;               ///< helper to define type of vector with element T
      };
  };

  struct BS_API_PLUGIN base_strategy_di : public base_strategy<double, double, int>
  {
  };
  struct BS_API_PLUGIN base_strategy_fi : public base_strategy<float, float, int>
  {
  };
  struct BS_API_PLUGIN base_strategy_mixi : public base_strategy<float, double, int>
  {
  };
  
#define BS_INST_STRAT(T)  \
  template class T<blue_sky::base_strategy_fi>;     \
  template class T<blue_sky::base_strategy_di>;     \
  template class T<blue_sky::base_strategy_mixi>;


#ifdef _MPI
  template <class fp_type, class i_type>
  struct BS_API_PLUGIN mpi_strategy
  {
    typedef fp_type                       item_t;             ///< short name for floating point type (old array_item_t)
    typedef i_type                        index_t;            ///< short name for index type
    typedef mpi_vector<item_t>            item_array_t;       ///< short name for array type (old array_t)
    typedef mpi_vector<index_t>           index_array_t;      ///< short name for type of array of indexes
    typedef matrix_base<item_array_t, index_array_t>
    matrix_t;           ///< short name for matrix type

    typedef mpi_csr_matrix <item_array_t, index_array_t>
    csr_matrix_t;       ///< short name for matrix type

    typedef std::vector<matrix_t>         matrix_array_t;     ///< short name for vector of matrices type
    typedef seq_vector<fp_type>           wksp_t;             ///< short name for workspace array

    template <typename T>
    struct vec
      {
        typedef mpi_vector <T>              type;               ///< helper to define type of vector with element T
      };

    typedef mpi_barrier                   barrier_t;          ///< short name for barrier
  };

  struct BS_API_PLUGIN mpi_strategy_di : public mpi_strategy <double, int>
  {
  };

  struct BS_API_PLUGIN mpi_strategy_fi : public mpi_strategy <float, int>
  {
  };
#endif
} // namespace blue_sky


#endif // #ifndef BS_LINEAR_SOLVERS_MATRIX_STRATEGIES_H_
