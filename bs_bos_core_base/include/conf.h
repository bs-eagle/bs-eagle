/** 
 * @file conf.h
 * @brief define default types  
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2011-02-22
 */
#ifndef CONF_7FBMZIFT

#define CONF_7FBMZIFT


//#ifdef BS_ARRAY_DEFAULT_TRAITS
//#undef BS_ARRAY_DEFAULT_TRAITS
//#define BS_ARRAY_DEFAULT_TRAITS  bs_nparray
//#endif //BS_ARRAY_DEFAULT_TRAITS

#include "bs_nparray.h"
#include "bs_array.h"

//! using this type for small data
#define t_int           int

//! this type should be using for indexing (matrix, arrays, ...)
#define t_long          long

//! this type should be using for single precision elements (matrix, fixed arrays) 
#define t_float         double

//! this type should be using for double precision elements (calculed arrays, arrays of unknowns) 
#define t_double        double

#define t_uint          unsigned long
#define t_ulong         unsigned long

#define v_int          bs_array<t_int, bs_nparray> 
#define v_long         bs_array<t_long, bs_nparray> 
#define v_uint         bs_array<t_uint, bs_nparray> 
#define v_ulong        bs_array<t_ulong, bs_nparray> 
#define v_float        bs_array<t_float, bs_nparray> 
#define v_double       bs_array<t_double, bs_nparray> 

#define stdv_int       std::vector<t_int> 
#define stdv_long      std::vector<t_long> 
#define stdv_uint      std::vector<t_uint> 
#define stdv_ulong     std::vector<t_ulong> 
#define stdv_float     std::vector<t_float> 
#define stdv_double    std::vector<t_double> 

#define spv_int          smart_ptr<v_int, true> 
#define spv_long         smart_ptr<v_long, true> 
#define spv_uint         smart_ptr<v_uint, true> 
#define spv_ulong        smart_ptr<v_ulong, true> 
#define spv_float        smart_ptr<v_float, true> 
#define spv_double       smart_ptr<v_double, true> 

#endif /* end of include guard: CONF_7FBMZIFT */
