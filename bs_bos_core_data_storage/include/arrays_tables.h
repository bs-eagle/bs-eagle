#ifndef _ARRAYS_TABLES_H
#define _ARRAYS_TABLES_H
#include "arrays.h"

extern BS_API_PLUGIN int i_pool_default_values[ARR_TOTAL_INT];
extern BS_API_PLUGIN double d_pool_default_values[ARR_TOTAL_DOUBLE];
extern BS_API_PLUGIN double c_pool_default_values[COMP_ARR_TOTAL_DOUBLE];
extern BS_API_PLUGIN int i_pool_sizes[ARR_TOTAL_INT * 6];
extern BS_API_PLUGIN int d_pool_sizes[ARR_TOTAL_DOUBLE * 6];
extern BS_API_PLUGIN const char *int_names_table [ARR_TOTAL_INT];
extern BS_API_PLUGIN const char *double_names_table [ARR_TOTAL_DOUBLE];
extern BS_API_PLUGIN const char *comp_names_table [COMP_ARR_TOTAL_DOUBLE];

#endif //_ARRAYS_TABLES_H
