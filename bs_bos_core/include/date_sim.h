#ifndef DATE_SIM_H
#define DATE_SIM_H
/*!
  \file date_sim.h
  \brief typedef of date structure
*/

//! Date structure
typedef struct tm DATE_SIM;

inline bool operator < (const tm &left, const tm &right)
{
  if (left.tm_year<right.tm_year) return true;
  if (left.tm_year>right.tm_year) return false;
  else
    {
      if (left.tm_mon<right.tm_mon) return true;
      if (left.tm_mon>right.tm_mon) return false;
      else
        {
          if (left.tm_mday<right.tm_mday) return true;
          if (left.tm_mday>right.tm_mday) return false;
        }
    }
  return false;
}

#endif // DATE_SIM_H
