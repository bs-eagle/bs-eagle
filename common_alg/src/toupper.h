#ifndef TOUPPER_N49LZ6PA

#define TOUPPER_N49LZ6PA


#include <ctype.h>

inline void 
locale_ucase (char *s)
{
  char *ptr = s;
  for (; *ptr != '\0'; ++ptr)
    *ptr = toupper (*ptr);
}

#endif /* end of include guard: TOUPPER_N49LZ6PA */
