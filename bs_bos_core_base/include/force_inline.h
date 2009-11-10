/**
 * \file force_inline.h
 * \brief cross-platform forceinline definition
 * \author Sergey Miryanov
 * \date 30.12.2008
 * */
#ifndef BS_FORCE_INLINE_H_
#define BS_FORCE_INLINE_H_

#if defined(_WIN32) && defined(_MSC_VER)
  #define BS_FORCE_INLINE __forceinline
#else // UNIX
  #define BS_FORCE_INLINE __attribute__ ((always_inline))
#endif // defined(_WIN32) && defined(_MSC_VER)

#endif  // #ifndef BS_FORCE_INLINE_H_
