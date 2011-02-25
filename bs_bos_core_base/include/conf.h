/** 
 * @file conf.h
 * @brief define default types  
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2011-02-22
 */
#ifndef CONF_7FBMZIFT

#define CONF_7FBMZIFT

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


#define v_int          bs_array<t_int>
#define v_long         bs_array<t_long>
#define v_uint         bs_array<t_uint>
#define v_ulong        bs_array<t_ulong>
#define v_float        bs_array<t_float>
#define v_double       bs_array<t_double>

#define spv_int          smart_ptr<v_int, true>
#define spv_long         smart_ptr<v_long, true>
#define spv_uint         smart_ptr<v_uint, true>
#define spv_ulong        smart_ptr<v_ulong, true>
#define spv_float        smart_ptr<v_float, true>
#define spv_double       smart_ptr<v_double, true>



#endif /* end of include guard: CONF_7FBMZIFT */
