#ifndef __MACROSES_H__
#define __MACROSES_H__
/*!
 * \file macroses.h
 * \brief usefull macroses
 * \author Borschuk Oleg
 * \date 2005-09-26
 */
#include "memory_macroses.h"
#include "debug_macroses.h"

//! calculate vector from begin and end points, N_S -- start point, N_E -- end point, V -- output vector
#define CALC_VECTOR(N_S,N_E,V)   (V)[0] = (N_E).array[0] - (N_S).array[0];(V)[2] = (N_E).array[2] - (N_S).array[2];(V)[1] = (N_E).array[1] - (N_S).array[1];

//! calculate triange area, V1 -- first vector, V2 -- second vector
#define CALC_TRIANGLE_AREA(V1,V2)  (0.5 * fabs((V1)[1] * (V2)[2] - (V1)[2] * (V2)[1] + (V1)[2] * (V2)[0] - (V1)[0] * (V2)[2] + (V1)[0] * (V2)[1] - (V1)[1] * (V2)[0]))

//! calculate vect mult
#define CALC_VECT_MULT(V1,V2,R)                         \
        R[0] = (V1)[1] * (V2)[2] - (V1)[2] * (V2)[1];   \
        R[1] = (V1)[2] * (V2)[0] - (V1)[0] * (V2)[2];   \
        R[2] = (V1)[0] * (V2)[1] - (V1)[1] * (V2)[0];

//! calculate vect norm
#define CALC_VECT_NORM(V) sqrt((V)[0] * (V)[0] + (V)[1] * (V)[1] + (V)[2] * (V)[2])
//! calculate geometric average
#define CALC_GEOM_AVERAGE(D1,D2)  ((D1) * (D2) / ((D1) + (D2)))

//! calculate cell ijk index from linear index IND -- linear index, NX, NY, NZ -- dimens, I, J, K -- output indexes
#define CALC_IJK_FROM_INDEX(IND,NX,NY,NZ,I,J,K)                         \
  (K) = (int)((IND) / ((NX) * (NY)));                                   \
  (IND) -= (K) * (NX) * (NY);                                           \
  (J) = (int)((IND) / (NX));                                            \
  (I) = (IND) - (J) * (NX);

//! calculate cell ijk index from linear index IND -- linear index, NX, NY, NZ -- dimens, I, J, K -- output indexes
//Defined by GataullinT
#define CALC_INDEX_FROM_IJK(IND,NX,NY,NZ,I,J,K)                         \
IND = I + J * NX + K * NX * NY                                          \
 

/////////////////// UFASOLVER VERSION AND DEMO ////////////////////////////////////
#if defined(_UFASOLVER_DEMO_VERSION_) || defined(_SGM_KERNEL_DEMO_)
//! max number of elements in mesh
#define DEMO_MAX_NUMBER_OF_ELEMENTS 1000
//! max number of wells in reservoir status
#define DEMO_MAX_NUMBER_OF_WELLS    10
//! check mesh elements
#define CHECK_DEMO_VERSION(N)                                           \
if ((N) > DEMO_MAX_NUMBER_OF_ELEMENTS)                                  \
    {                                                                   \
      rep->print (LOG_INIT_SECTION, LOG_CRIT,                           \
                  GET_TEXT ("\nSorry, but this DEMO version not support meshes larger than %d grid blocks.\n"), DEMO_MAX_NUMBER_OF_ELEMENTS);              \
      rep->print (LOG_INIT_SECTION, LOG_CRIT,                           \
                  GET_TEXT ("For purchase full version please email on bos@ngt.ru \nor see http://www.ngt.ru/soft/bos.html for more information.\n"));   \
      return -13;                                                       \
    }
//! check number of wells
#define CHECK_DEMO_WELLS(OBJECT,RET)                                    \
  if ((OBJECT)->demo_well_count > DEMO_MAX_NUMBER_OF_WELLS)             \
    {                                                                   \
      rep->print (LOG_INIT_SECTION, LOG_CRIT,                           \
                  GET_TEXT ("\nSorry, but this DEMO version not support number of wells greater than %d.\n"), DEMO_MAX_NUMBER_OF_WELLS);              \
      rep->print (LOG_INIT_SECTION, LOG_CRIT,                           \
                  GET_TEXT ("For purchase full version please email on bos@ngt.ru \nor see http://www.ngt.ru/soft/bos.html for more information.\n"));   \
      return (RET);                                                     \
    }
#endif // _UFASOLVER_DEMO_VERSION_

//! Minimum of a and b
#define MIN(a,b) ((a)<(b)?(a):(b))
//! Maximum of a and b
#define MAX(a,b) ((a)<(b)?(b):(a))

#endif //__MACROSES_H__
