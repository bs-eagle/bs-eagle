/**
 *       \file  date_sim.h
 *      \brief  Typedef of date structure
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef DATE_SIM_H
#define DATE_SIM_H

//! Date structure
typedef struct tm DATE_SIM;

/**
 * \brief  Compares two dates
 * \param  left
 * \param  right
 * \return True if right > left
 * */
inline bool 
operator < (const tm &left, const tm &right)
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
