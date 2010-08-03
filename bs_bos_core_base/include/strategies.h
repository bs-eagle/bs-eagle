/**
* \file   call_helper.cpp
* \brief  Strategies for linear solvers
* \author Miryanov Sergey
* \date 2008-04-03
*/

#ifndef BS_LINEAR_SOLVERS_MATRIX_STRATEGIES_H_
#define BS_LINEAR_SOLVERS_MATRIX_STRATEGIES_H_

#include <string>

#ifdef _MPI
#include <mpi.h>
#include "mpi_vector.h"
#endif

#include "shared_vector.h"
#include "strategy_iface.h"

namespace blue_sky
{
  //template <typename array_item_t, typename array_index_t>
  //template <typename fp_vector_type, typename i_vector_type>
  //class matrix_iface;

  //template <typename fp_vector_type, typename i_vector_type, typename fp_storage_vector_type>
  //class bcsr_matrix_iface;

  //template <typename array_item_t, typename array_index_t>
  //class BS_API_PLUGIN bcsr_matrix;

  //template <typename array_item_t, typename array_index_t>
  //class BS_API_PLUGIN mpi_csr_matrix;

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
#if 0
  template <class rhs_fp_type, class fp_type, class i_type>
  struct BS_API_PLUGIN base_strategy
  {
    typedef fp_type                                           item_t;             ///< short name for floating point type (old array_item_t)
    typedef i_type                                            index_t;            ///< short name for index type
    typedef seq_vector<item_t>                                item_array_t;       ///< short name for array type (old array_t)
    typedef seq_vector<index_t>                               index_array_t;      ///< short name for type of array of indexes
    typedef rhs_fp_type		                                    rhs_item_t;         ///< short name for type of rhs
    typedef seq_vector<rhs_item_t>	                          rhs_item_array_t;   ///< short name for rhs array type

    typedef matrix_base<rhs_item_array_t, index_array_t>      matrix_t;           ///< short name for matrix type
    typedef bcsr_matrix<rhs_item_array_t, index_array_t>      csr_matrix_t;       ///< short name for csr_matrix type (it's fail way today)

    typedef std::vector<matrix_t>                             matrix_array_t;     ///< short name for vector of matrices type
    typedef seq_vector<fp_type>                               wksp_t;             ///< short name for workspace array

    template <typename T>
    struct vec
      {
        typedef seq_vector <T>              type;               ///< helper to define type of vector with element T
      };

    typedef seq_barrier                   barrier_t;          ///< short name for barrier
  };
#endif //0


  template <class fp_type, class i_type, class fp_storage_type>
  struct BS_API_PLUGIN base_strategy : public strategy_iface
  {
    //! type for storing fp vectors (not in matrix)
    typedef shared_vector<fp_type>                 fp_vector_type;
    //! type for storing i vectors 
    typedef shared_vector<i_type>                  i_vector_type;
    //! type for storing fp vectors in matrices
    typedef shared_vector<fp_storage_type>         fp_storage_vector_type;

    //typedef matrix_iface<fp_vector_type, i_vector_type>      matrix_t;           ///< short name for matrix type

    typedef fp_type             fp_type_t;
    typedef i_type              i_type_t;
    typedef fp_storage_type     fp_storage_type_t;
    //typedef std::vector<matrix_t>                             matrix_array_t;     ///< short name for vector of matrices type
    
    //typedef seq_vector<fp_type>                               wksp_t;             ///< short name for workspace array

    template <typename T>
    struct vec
      {
        typedef shared_vector <T>              type;               ///< helper to define type of vector with element T
      };

    typedef seq_barrier                   barrier_t;          ///< short name for barrier

  };


  struct BS_API_PLUGIN base_strategy_did : public base_strategy<double, int, double>
  {
    public:
      static std::string get_suf () {return std::string ("did");}
  };
  struct BS_API_PLUGIN base_strategy_fif : public base_strategy<float, int, float>
  {
    public:
      static std::string get_suf () {return std::string ("fif");}
  };
  struct BS_API_PLUGIN base_strategy_dif : public base_strategy<double, int, float>
  {
    public:
      static std::string get_suf () {return std::string ("dif");}
  };

  struct BS_API_PLUGIN base_strategy_dld : public base_strategy<double, long, double>
  {
    public:
      static std::string get_suf () {return std::string ("dld");}
  };
  struct BS_API_PLUGIN base_strategy_flf : public base_strategy<float, long, float>
  {
    public:
      static std::string get_suf () {return std::string ("flf");}
  };
  struct BS_API_PLUGIN base_strategy_dlf : public base_strategy<double, long, float>
  {
    public:
      static std::string get_suf () {return std::string ("dlf");}
  };

  
#define BS_INST_STRAT(T)  \
  template class T<blue_sky::base_strategy_fif>;        \
  template class T<blue_sky::base_strategy_did>;        \
  template class T<blue_sky::base_strategy_dif>;        \
  template class T<blue_sky::base_strategy_flf>;        \
  template class T<blue_sky::base_strategy_dld>;        \
  template class T<blue_sky::base_strategy_dlf>;


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
    typedef shared_vector<fp_type>           wksp_t;             ///< short name for workspace array

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
