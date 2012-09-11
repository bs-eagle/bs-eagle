/** 
 * @file conf.h
 * @brief define default types  
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2011-02-22
 */
#ifndef CONF_7FBMZIFT
#define CONF_7FBMZIFT

#ifdef BS_ARRAY_DEFAULT_TRAITS
#undef BS_ARRAY_DEFAULT_TRAITS
#endif
#define BS_ARRAY_DEFAULT_TRAITS bs_nparray

#include "bs_nparray.h"
#include "bs_npvec.h"
#include "bs_npvec_shared.h"
#include "bs_array.h"

//! using this type for small data
typedef int           t_int;

//! this type should be using for indexing (matrix, arrays, ...)
typedef long long     t_long;

#ifndef T_FLOAT_IS_DOUBLE
#define T_FLOAT_IS_DOUBLE 1
#endif

#ifdef T_FLOAT_IS_DOUBLE
//! this type should be using for single precision elements (matrix, fixed arrays) 
typedef double        t_float;
#else
typedef float         t_float;
#endif

//! this type should be using for double precision elements (calculed arrays, arrays of unknowns) 
typedef double        t_double;

typedef unsigned long               t_uint;
typedef unsigned long               t_ulong;

#define BS_EAGLE_ARRAY_TRAITS blue_sky::bs_npvec_shared
typedef blue_sky::bs_array< t_int    , BS_EAGLE_ARRAY_TRAITS >    v_int;
typedef blue_sky::bs_array< t_long   , BS_EAGLE_ARRAY_TRAITS >   v_long;
typedef blue_sky::bs_array< t_uint   , BS_EAGLE_ARRAY_TRAITS >   v_uint;
typedef blue_sky::bs_array< t_ulong  , BS_EAGLE_ARRAY_TRAITS >  v_ulong;
typedef blue_sky::bs_array< t_float  , BS_EAGLE_ARRAY_TRAITS >  v_float;
typedef blue_sky::bs_array< t_double , BS_EAGLE_ARRAY_TRAITS > v_double;

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
namespace blue_sky {

std::string get_long_type ();
namespace python {

void export_type_helper();

} // namespace python
} // namespace blue_sky



#endif /* end of include guard: CONF_7FBMZIFT */
