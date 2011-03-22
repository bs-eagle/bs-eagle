/** 
 * @file conf.h
 * @brief define default types  
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2011-02-22
 */
#ifndef CONF_7FBMZIFT
#define CONF_7FBMZIFT

#include "bs_nparray.h"
#include "bs_array.h"

namespace blue_sky {

template< class T >
struct bs_array_default_traits< T, bs_array_default_traits_tag > : public bs_array_base2traits< T, bs_nparray > {};

}

//! using this type for small data
typedef int           t_int;

//! this type should be using for indexing (matrix, arrays, ...)
typedef long          t_long;

//! this type should be using for single precision elements (matrix, fixed arrays) 
typedef double        t_float;

//! this type should be using for double precision elements (calculed arrays, arrays of unknowns) 
typedef double        t_double;

typedef unsigned long               t_uint;
typedef unsigned long               t_ulong;

typedef blue_sky::bs_array< t_int >           v_int;
typedef blue_sky::bs_array< t_long >          v_long;
typedef blue_sky::bs_array< t_uint >          v_uint;
typedef blue_sky::bs_array< t_ulong >         v_ulong;
typedef blue_sky::bs_array< t_float >         v_float;
typedef blue_sky::bs_array< t_double >        v_double;

typedef std::vector< t_int >        stdv_int;
typedef std::vector< t_long >       stdv_long;
typedef std::vector< t_uint >       stdv_uint;
typedef std::vector< t_ulong >      stdv_ulong;
typedef std::vector< t_float >      stdv_float;
typedef std::vector< t_double >     stdv_double;

typedef blue_sky::smart_ptr< v_int >          spv_int;
typedef blue_sky::smart_ptr< v_long >         spv_long;
typedef blue_sky::smart_ptr< v_uint >         spv_uint;
typedef blue_sky::smart_ptr< v_ulong >        spv_ulong;
typedef blue_sky::smart_ptr< v_float >        spv_float;
typedef blue_sky::smart_ptr< v_double >       spv_double;

#endif /* end of include guard: CONF_7FBMZIFT */
