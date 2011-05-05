#ifndef __MEMORY_MACROSES_H_
#define __MEMORY_MACROSES_H_
/*!
 * \file memory_macroses.h
 * \brief usefull macroses for use in solver
 * \author Borschuk Oleg
 * \date 2006-06-19
 */

//! fill array A -- array, N1 -- start, N2 -- end, V -- value
#define FI_FILL_ARRAY(A,N1,N2,V)                                        \
  for (t_long iii = (N1); iii < (N2); ++iii)                            \
    {                                                                   \
      (A)[iii] = (V);                                                   \
    }
//! decrement array by value A -- array, N1 -- start, N2 -- end, V -- value
#define FI_DECR_ARRAY(A,N1,N2,V)                                        \
  for (t_long iii = (N1); iii < (N2); ++iii)                            \
    {                                                                   \
      (A)[iii] -= (V);                                                  \
    }

#endif //__MEMORY_MACROSES_H_
