#ifndef __DEBUG_MACROSES_H_
#define __DEBUG_MACROSES_H_
/*!
 * \file debug_macroses.h
 * \brief usefull macroses for debug
 * \author Borschuk Oleg
 * \date 2006-06-19
 */

#include <stdlib.h>
#include <stdio.h>
//#include "sysinfo.h"

///////////////////////////////////////////////////////////////////////////////////
//! author KarimovRF,
//! date 2005-10-06
//! check array A access for index I conforming with array size SZ.
//! If the given index is illegal (is less than zero or greater than the array size)
//! then the error_messager (see sysinfo.h) function is called.
//! Otherwise macro returns I-th element of A array
#ifdef MEM_DEBUG
#define CHA(A,I,SZ) ((((I) < 0) || ((I) > (SZ - 1))) ? (A)[error_messager (__FILE__, __LINE__, " Array bounds overrun!" )] : (A)[I])
#else
#define CHA(A,I,SZ) ((A)[I])
#endif

///////////////////////////////////////////////////////////////////////////////////
//! author SharipovTR,
//! date 2006-01-06
//! check array A access for index I conforming with array size SZ.
//! If the given index is illegal (is less than zero or greater than the array size)
//! then the error_messager (see sysinfo.h) function is called.
//! Otherwise macro returns I-th element of A array
#define ERROR_MSG() error_messager (__FILE__, __LINE__, " Array bounds overrun!" )

#define MEMORY_ERROR_REPORT_PRINT()                                           \
  {                                                                            \
    rep->print (                                                               \
      LOG_READ_SECTION, LOG_CRIT, "%s: %s\n",                                  \
      GET_TEXT ("Fatal error"),                                                \
      GET_TEXT ("not enough memory!"));                                        \
  }

#ifdef MEM_DEBUG
#define ITEM(A,I,J,SZI,SZJ) ((((I) < 0) || ((I) > (SZI - 1)) || ((J) < 0) || ((J) > (SZJ - 1))) ?    \
                            (A)[ERROR_MSG()][ERROR_MSG()] : (A)[I][J])
#else
#define ITEM(A,I,J,SZI,SZJ) ((A)[I][J])
#endif

//! author KarimovRF,
//! date 2005-10-06
//! checking specified pointer for being initialized. If not, error message
//! is generated, and ERR value is returned
#ifdef MEM_DEBUG
#define CH_PTR(PT)                                                          \
   if (!(PT))                                                                   \
     {                                                                          \
       error_messager (__FILE__, __LINE__, " null pointer!" );                  \
     }
#else //MEM_DEBUG
#define CH_PTR(PT)   ;
#endif //MEM_DEBUG    


/*!
 * \brief Base class for ASSERT macros implementation
 */
class asserter
  {
  public:
    asserter(int holds, const char *msg) : holds_(holds), msg_(msg) {}

    virtual ~asserter() {}

    virtual int handle(const char * /*file*/, int /*line*/) const
      {
        if (holds_)
          return true;

        enum {n = 2048};
        //static char buf[n];

//!	      sprintf(buf, "\n\tSource File: %s, Line: %d\n\tMessage:     %s\n\n",
        //   file, line, msg_ ? msg_ : "");

        //crit err

//!		  //print log fatal error
        //fatal_error_handler(buf);

        exit(-1);
      }

  protected:
    const int holds_;
    const char *msg_;

    asserter &operator = (const asserter &rhs);
  };

/**
 * BUG: Expansion of __LINE__ Macro Is Incorrect When Using /ZI
 * RESOLUTION: To work around this problem, use the Program Database option (/Zi)
 * instead of Edit and Continue (/ZI).
 *
 * Using an trick for avoid duplicate identifiers error in case:
 *      ASSERT(exp1); ASSERT(exp2);
 */
#define BOS_ASSERT(x)\
    {\
        asserter info(x, "Assertion failed: "#x);\
        info.handle(__FILE__, __LINE__);\
    }

#endif //__DEBUG_MACROSES_H_

