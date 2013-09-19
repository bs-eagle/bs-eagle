#ifndef PURE_MESH_H
#define PURE_MESH_H
/*!
	\file pure_mesh.h
	\brief This file declares definitions for pure meshe without bs
	\author Mark Khait
	\date 2012-03-22
 */

#include "csr_matrix.h"
#include "memory_macroses.h"
#include <stdio.h>
//#include <boost/array.hpp>
#include <vector>

#define BS_API_PLUGIN 
#define BS_ASSERT
#define bs_throw_exception printf
#define write_time_to_log printf
#define init_time
#ifndef CH_PTR
#   define CH_PTR
#endif // CH_PTR
#define BS_SP(T) T


//! using this type for small data
typedef int           t_int;

//! this type should be using for indexing (matrix, arrays, ...)
typedef long          t_long;
typedef double        t_double;
typedef float        t_float;

typedef int*          spv_int;
typedef long*         spv_long;
typedef float*       spv_float;
typedef double*       spv_double;


typedef std::vector< t_int >        stdv_int;
typedef std::vector< t_long >       stdv_long;
typedef std::vector< t_float >      stdv_float;
typedef std::vector< t_double >     stdv_double;

struct sp_hdm
  {
    t_float minpv;
    t_float minsv;
    t_float max_thickness;
    t_float darcy_constant;
    t_long n_elements;
    t_int nx, ny, nz;

    spv_int actnum_array;
    spv_float poro_array;
    spv_float ntg_array;
    spv_float multpv_array;
    spv_float permx_array;
    spv_float permy_array;
    spv_float permz_array;
    spv_float multx_array;
    spv_float multy_array;
    spv_float multz_array;

    spv_float zcorn_array;
    spv_float coord_array;
  };

struct sp_flux_conn_iface
  {
    csr_matrix conn_trans;
    spv_long matrix_block_idx_plus;
    spv_long matrix_block_idx_minus;
  };

typedef sp_hdm*             sp_hdm_t;
typedef sp_flux_conn_iface* sp_flux_conn_iface_t;

typedef void* sp_well_pool_t;
typedef void * table_iface;



    






#endif //PURE_MESH_H
