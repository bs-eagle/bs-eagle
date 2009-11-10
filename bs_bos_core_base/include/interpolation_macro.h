#ifndef __INTERPOLATION_MACRO_H
#define __INTERPOLATION_MACRO_H

/*!
  \file interpolation.h
  \brief Interpolation algotithms

  Interpolation algotithms realized as macros to reduce number of calls

*/

//! \brief Locate X0 in X (nreg x ntab) FORTRAN array by binary search algorithm
//! \param  X - (nreg x ntab) FORTRAN array
//! \param X0 - value to look up
//! \param N - array length
//! \param NREG - array dimension (so, it is NREG x N array)
//! \param REGNUM - region number
//! \param IL - local variable
//! \param IU - local variable
//! \param IM - local variable
//! \return IL - answer
#define BINARY_SEARCH_REG(X,X0,NREG,N,REGNUM,IU,IM,IL)     \
  do                                                    \
    {                                                   \
      (IL) = 0;                                         \
      (IU) = (N);                                       \
      while ((IL) != (IU))                              \
        {                                               \
          (IM) = ((IL) + (IU)) / 2;                     \
          if ((X)[(REGNUM) + (IM) * (NREG)] < (X0))     \
            (IL) = (IM) + 1;                            \
          else                                          \
            (IU) = (IM);                                \
        }                                               \
    }                                                   \
  while (0)

// New interpolation macros for data array
#define BINARY_SEARCH(X,X0,N,IU,IM,IL)                  \
  do                                                    \
    {                                                   \
      (IL) = 0;                                         \
      (IU) = (N);                                       \
      while ((IL) != (IU))                              \
        {                                               \
          (IM) = ((IL) + (IU)) / 2;                     \
          if ((X)[(IM)] < (X0))                         \
            (IL) = (IM) + 1;                            \
          else                                          \
            (IU) = (IM);                                \
        }                                               \
    }                                                   \
  while (0)

// Interpolation macros used in capp_equil.cpp
#define UP_BINARY_SEARCH(X,X0,N,IU,IM,IL)               \
  do                                                    \
    {                                                   \
      (IL) = 0;                                         \
      (IU) = (N);                                       \
      while ((IL) != (IU))                              \
        {                                               \
          (IM) = ((IL) + (IU)) / 2;                     \
          if ((X)[(IM)] > (X0))                         \
            (IL) = (IM) + 1;                            \
          else                                          \
            (IU) = (IM);                                \
        }                                               \
    }                                                   \
  while (0)


//! \brief Compute interval for linear interpolation
//! \param IL - computed interval on segment [0,N]
//! \param NREG - array dimension
//! \param N - lenght of interval
//! \param REGNUM - region number
//! \param I1 - result (bounds for linear interpolation)
//! \param I2 - result (bounds for linear interpolation)
//! \return I1, I2
#define COMPUTE_INTERPOLATION_INTERVAL_REG(IL,NREG,N,REGNUM,I1,I2) \
((IL) == (N)                                                    \
 ?    ((I1) = (REGNUM) + ((IL) - 2) * (NREG),                   \
       (I2) = (REGNUM) + ((IL) - 1) * (NREG))                   \
 : ((IL) == 0                                                   \
    ? ((I1) = (REGNUM) + ((IL)    ) * (NREG),                   \
       (I2) = (REGNUM) + ((IL) + 1) * (NREG))                   \
    : ((I1) = (REGNUM) + ((IL) - 1) * (NREG),                   \
       (I2) = (REGNUM) + ((IL)    ) * (NREG))))


// New interpolation macros for data array
#define COMPUTE_INTERPOLATION_INTERVAL(IL,N,I1,I2) \
((IL) == (N)                    \
 ?    ((I1) = ((IL) - 2),       \
       (I2) = ((IL) - 1))       \
 : ((IL) == 0                   \
    ? ((I1) = ((IL)    ),       \
       (I2) = ((IL) + 1))       \
    : ((I1) = ((IL) - 1),       \
       (I2) = ((IL)    ))))

//! \brief Compute picewise linear interpolation (with constant extrapolation)\n
//! of function specified by table of values Y (nreg x ntab FORTRAN array)\n
//! in points X (nreg x ntab FORTRAN array)
//! \param  X - (nreg x ntab) FORTRAN array
//! \param  Y - (nreg x ntab) FORTRAN array
//! \param X0 - value to look up
//! \param Y0 - answer
//! \param N - array length
//! \param NREG - array dimension (so, it is NREG x N array)
//! \param REGNUM - region number
//! \return Y0 - answer
#define COMPUTE_LINEAR_CONSTANT_INTERPOLATION_REG(X,Y,X0,Y0,NREG,N,REGNUM) \
  do                                                                    \
    {                                                                   \
      int il, iu, im;                                                   \
      /* Locate X0 in array X by binary search algorithm */             \
      BINARY_SEARCH_REG (X, X0, NREG, N, REGNUM, iu, im, il);              \
                                                                        \
      if (il == (N))                                                    \
        /* Constant extrapolation */                                    \
        (Y0) = (Y)[(REGNUM) + ((N) - 1) * (NREG)];                      \
      else if (il == 0)                                                 \
        /* Constant extrapolation */                                    \
        (Y0) = (Y)[(REGNUM)];                                           \
      else                                                              \
        {                                                               \
          /* Linear interpolation */                                    \
          im = (REGNUM) + (il - 1) * (NREG);                            \
          iu = (REGNUM) + il * (NREG);                                  \
          (Y0) = (Y)[im]                                                \
                 + ((X0) - (X)[im]) * ((Y)[iu] - (Y)[im])               \
                                    / ((X)[iu] - (X)[im]);              \
        }                                                               \
    }                                                                   \
  while (0)

// interpolate internal interval lineary and extrapolate with fixet derivates
#define INTERNAL_INTERPOLATION_AND_FIX_D_EXTRAPOLATION(X,Y,X0,Y0,N)             \
  do{                                                                           \
      if (il == (N))                                                            \
        {                                                                       \
          if ((N) == 1)                                                         \
            /* Constant extrapolation */                                        \
            (Y0) = (Y)[((N) - 1)];                                              \
          else                                                                  \
            (Y0) = (Y)[(N) - 2]                                                 \
                 + ((X0) - (X)[(N) - 2]) * ((Y)[(N) - 1] - (Y)[(N) - 2])        \
                                         / ((X)[(N) - 1] - (X)[(N) - 2]);       \
        }                                                                       \
      else if (il == 0)                                                         \
        {                                                                       \
          if ((N) == 1)                                                         \
            /* Constant extrapolation */                                        \
            (Y0) = (Y)[0];                                                      \
          else                                                                  \
            (Y0) = (Y)[0]                                                       \
                 + ((X0) - (X)[0]) * ((Y)[1] - (Y)[0])                          \
                                         / ((X)[1] - (X)[0]);                   \
        }                                                                       \
      else                                                                      \
        {                                                                       \
          /* Linear interpolation */                                            \
          im = (il - 1);                                                        \
          iu = il;                                                              \
          (Y0) = (Y)[im]                                                        \
                 + ((X0) - (X)[im]) * ((Y)[iu] - (Y)[im])                       \
                                    / ((X)[iu] - (X)[im]);                      \
        }                                                                       \
   }while (0)

// interpolate internal interval lineary and extrapolate with fixet derivates
#define INTERNAL_INTERPOLATION_AND_FIX_D_EXTRAPOLATION_DERIV(X,Y,X0,DY0,N)      \
  do{                                                                           \
      if (il == (N))                                                            \
        {                                                                       \
          if ((N) == 1)                                                         \
            /* Constant extrapolation */                                        \
            (DY0) = 0;                                                          \
          else                                                                  \
            (DY0) = ((Y)[(N) - 1] - (Y)[(N) - 2])                               \
                                         / ((X)[(N) - 1] - (X)[(N) - 2]);       \
        }                                                                       \
      else if (il == 0)                                                         \
        {                                                                       \
          if ((N) == 1)                                                         \
            /* Constant extrapolation */                                        \
            (DY0) = 0;                                                          \
          else                                                                  \
            (DY0) = ((Y)[1] - (Y)[0]) / ((X)[1] - (X)[0]);                      \
        }                                                                       \
      else                                                                      \
        {                                                                       \
          /* Linear interpolation */                                            \
          im = (il - 1);                                                        \
          iu = il;                                                              \
          (DY0) = ((Y)[iu] - (Y)[im]) / ((X)[iu] - (X)[im]);                    \
        }                                                                       \
   }while (0)


#define INTERNAL_INTERPOLATION(X,Y,X0,Y0,N)                             \
   do {                                                                 \
      if (il == (N))                                                    \
        /* Constant extrapolation */                                    \
        (Y0) = (Y)[((N) - 1)];                                          \
      else if (il == 0)                                                 \
        /* Constant extrapolation */                                    \
        (Y0) = (Y)[0];                                                  \
      else                                                              \
        {                                                               \
          /* Linear interpolation */                                    \
          im = (il - 1);                                                \
          iu = il;                                                      \
          (Y0) = (Y)[im]                                                \
                 + ((X0) - (X)[im]) * ((Y)[iu] - (Y)[im])               \
                                    / ((X)[iu] - (X)[im]);              \
        }                                                               \
   }while (0)

// New interpolation macros for data array
#define COMPUTE_LINEAR_CONSTANT_INTERPOLATION(X,Y,X0,Y0,N)              \
  do                                                                    \
    {                                                                   \
      int il, iu, im;                                                   \
      /* Locate X0 in array X by binary search algorithm */             \
      BINARY_SEARCH (X, X0, N, iu, im, il);                             \
                                                                        \
      INTERNAL_INTERPOLATION (X, Y, X0, Y0, N);                         \
    }                                                                   \
  while (0)

// New interpolation macros for data array
#define COMPUTE_INTERPOLATION_AND_DERIVATIVES(X,Y,DY,X0,Y0,DY0,N)       \
  do                                                                    \
    {                                                                   \
      int il, iu, im;                                                   \
      /* Locate X0 in array X by binary search algorithm */             \
      BINARY_SEARCH (X, X0, N, iu, im, il);                             \
      /* function interpolation */                                      \
      INTERNAL_INTERPOLATION (X, Y, X0, Y0, N);                         \
      /* derivatives interpolation */                                   \
      INTERNAL_INTERPOLATION (X, DY, X0, DY0, N - 1);                   \
    }                                                                   \
  while (0)


#define INTERNAL_CONST_INTERPOLATION(X,Y,X0,Y0,N)                       \
   do {                                                                 \
      if (il == (N))                                                    \
        /* Constant extrapolation */                                    \
        (Y0) = 0;                                                       \
      else if (il == 0)                                                 \
        /* Constant extrapolation */                                    \
        (Y0) = 0;                                                       \
      else                                                              \
        {                                                               \
          /* Linear interpolation */                                    \
          im = (il - 1);                                                \
          iu = il;                                                      \
          (Y0) = ((Y)[iu] - (Y)[im]) / ((X)[iu] - (X)[im]);             \
        }                                                               \
  }while (0)
// compute linear interpolation and constant derivation
// X --- X array
// Y --- Y array
#define COMPUTE_LIN_INTERP_AND_CONST_DERIV(X,Y,X0,Y0,DY0,N)             \
  do                                                                    \
    {                                                                   \
      int il, iu, im;                                                   \
      /* Locate X0 in array X by binary search algorithm */             \
      BINARY_SEARCH (X, X0, N, iu, im, il);                             \
      /* function interpolation */                                      \
      INTERNAL_INTERPOLATION (X, Y, X0, Y0, N);                         \
      /* derivatives interpolation */                                   \
      INTERNAL_CONST_INTERPOLATION (X, Y, X0, DY0, N);                  \
    }                                                                   \
  while (0)

// compute interpolation and constant derivation for phase permeabilities
// X --- X array
// Y --- Y array
//
#if 0
#define COMPUTE_INTERP_PHASE_PERM(X,Y,X0,Y0,DY0,N)                      \
  do                                                                    \
    {                                                                   \
      int izero;                                                        \
      for (izero = 0; (Y)[izero + 1] == 0; izero++)                     \
        ;                                                               \
      /* function interpolation */                                      \
      if (izero == (il - 1))                                            \
        (Y0) = (Y)[izero + 1]                                           \
                  * ((X0) - (X)[izero]) * ((X0) - (X)[izero])           \
                  / ((X)[izero + 1] - (X)[izero])                       \
                  / ((X)[izero + 1] - (X)[izero]);                      \
      else                                                              \
        INTERNAL_INTERPOLATION (X, Y, X0, Y0, N);                       \
      /* derivatives interpolation */                                   \
      if ((izero == (il - 1)) && (X0) != (X)[izero + 1])                \
        (DY0) = 2 * (Y)[izero + 1]                                      \
                       * ((X0) - (X)[izero])                            \
                       / ((X)[izero + 1] - (X)[izero])                  \
                       / ((X)[izero + 1] - (X)[izero]);                 \
      else if (il == (N))                                               \
        (DY0) = ((Y)[(N) - 1] - (Y)[(N) - 2])                           \
              / ((X)[(N) - 1] - (X)[(N) - 2]);                          \
      else                                                              \
        INTERNAL_CONST_INTERPOLATION (X, Y, X0, DY0, N);                \
    }                                                                   \
  while (0)
#else
#define COMPUTE_INTERP_PHASE_PERM(X,Y,X0,Y0,DY0,N)                      \
  do                                                                    \
    {                                                                   \
      INTERNAL_INTERPOLATION (X, Y, X0, Y0, N);                       \
      INTERNAL_CONST_INTERPOLATION (X, Y, X0, DY0, N);                \
    }                                                                   \
  while (0)
#endif
// compute interpolation and constant derivation for oil permeabilities
// X --- X array
// Y --- Y array

#define COMPUTE_INTERP_OIL_PERM(X,Y,X0,Y0,DY0,N)                        \
  do                                                                    \
    {                                                                   \
      int izero;                                                        \
      for (izero = (N) - 1; (Y)[izero - 1] == 0; izero--)               \
        ;                                                               \
      /* function interpolation */                                      \
      if (izero == il)                                                  \
        (Y0) = (Y)[izero - 1]                                           \
                  * ((X0) - (X)[izero]) * ((X0) - (X)[izero])           \
                  / ((X)[izero - 1] - (X)[izero])                       \
                  / ((X)[izero - 1] - (X)[izero]);                      \
      else                                                              \
        INTERNAL_INTERPOLATION (X, Y, X0, Y0, N);                       \
      /* derivatives interpolation */                                   \
      if ((izero == il) && (X0) != (X)[izero - 1])                      \
        (DY0) = 2 * (Y)[izero - 1]                                      \
                       * ((X0) - (X)[izero])                            \
                       / ((X)[izero - 1] - (X)[izero])                  \
                       / ((X)[izero - 1] - (X)[izero]);                 \
      else if (il == 0)                                                 \
        (DY0) = ((Y)[1] - (Y)[0])                                       \
              / ((X)[1] - (X)[0]);                                      \
      else                                                              \
        INTERNAL_CONST_INTERPOLATION (X, Y, X0, DY0, N);                \
    }                                                                   \
  while (0)

#define COMPUTE_LIN_INTERP(X1,Y1,X2,Y2,VX)  ((Y1) + ((VX) - (X1)) * ((Y2) - (Y1)) / ((X2) - (X1)))


#endif // __INTERPOLATION_MACRO_H
