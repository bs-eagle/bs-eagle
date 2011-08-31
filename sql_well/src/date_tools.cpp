/*!
  \file date_tools.cpp
  \brief date_tools - date handling tools
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "date_tools.h"
#include "main_def.h"
#include "localization.h"
//#include "error_num_def.h"
#include "timer.h"


/*! \brief Nonzero if YEAR is a leap year (every 4 years,\n
except every 100th isn't, and every 400th is).  */
#define IS_LEAP_YEAR(YEAR) \
  ((YEAR) % 4 == 0 && ((YEAR) % 100 != 0 || (YEAR) % 400 == 0))

 
/*!
  \brief Date in format %d %s %d (1 jan 2000) (reads from text buffer BUF\n
         into structure START
  \return -1 - in case of error\n
           0 - in case of success.
*/
int
key_read_date_ecl (const char *buf, int &day, int &month, int &year)
{
  char string[512];
  char *ptr_start = 0;
  char *ptr_end = 0;
  int i;
  
  if (!buf)
    return -3;
  
  // Need to add reading check
  if (sscanf (buf, "%d %s %d", &day, string, &year) != 3) 
    return -1;
    
  locale_ucase(string);
  ptr_start = string;
  
  while (*ptr_start == '\'' || *ptr_start == '\"')
    ++ptr_start;
  for (i = 0; ptr_start[i] != '\0'; ++i)
    ;
  ptr_end = &(ptr_start[i]) - 1;
  while (*ptr_end == '\'' || *ptr_end == '\"')
    --ptr_end;
  *(ptr_end + 1) = '\0';    
  if (!strcmp (ptr_start, "JAN"))
    month = 1;
  else if (!strcmp (ptr_start, "FEB"))
    month = 2;
  else if (!strcmp (ptr_start, "MAR"))
    month = 3;
  else if (!strcmp (ptr_start, "APR"))
    month = 4;
  else if (!strcmp (ptr_start, "MAY"))
    month = 5;
  else if (!strcmp (ptr_start, "JUN"))
    month = 6;
  else if (!strcmp (ptr_start, "JLY") || !strcmp (ptr_start, "JUL"))
    month = 7;
  else if (!strcmp (ptr_start, "AUG"))
    month = 8;
  else if (!strcmp (ptr_start, "SEP"))
    month = 9;
  else if (!strcmp (ptr_start, "OCT"))
    month = 10;
  else if (!strcmp (ptr_start, "NOV"))
    month = 11;
  else if (!strcmp (ptr_start, "DEC"))
    month = 12;
  else
    return -4;
  return 0;
}

int 
key_read_time (const char *buf, date_sim &tm)
{
  int h, m;
  double s;
  
  if (sscanf (buf, "%d:%d:%lf", &h, &m, &s) != 3)
    {
      fprintf (stderr, "Error: invalid time format: %s\n", buf);
      return -1;
    }
 
  // Check month and year
  
  if (h < 0 || h > 23)
    {
      fprintf (stderr, "Error: hours should be 0..23: %s\n", buf);
      return -1;
    }
  
  if (m < 0 || m > 59)
    {
      fprintf (stderr, "Error: minutes should be 0..59: %s\n", buf);
      return -1;
    } 
  if (s < 0 || s >= 60.0)
    {
      fprintf (stderr, "Error: seconds should be 0..59: %s\n", buf);
      return -1;
    } 
  tm = ((date_sim)h + ((date_sim)m + (date_sim)s / 60.0) / 60.0) / 24.0;
  return 0;
}

double
ymd2d (int year, int month, int day)
{
  int days_in_current_month;
  static const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int i;
  double date;
  
  // Check month and year
  
  if (year < 1900)
    {
      fprintf (stderr, "Error: years earlier 1900 not supported\n");
      return -1;
    }
  
  if (month < 1 || month > 12)
    {
      fprintf (stderr, "Error: invalid (nonexistent) date\n");
      return -1;
    } 
  
  date = 0;
  
  for (i = 1900; i < year; i++)
    {
    if (IS_LEAP_YEAR (i))
      date += 366;
    else
      date += 365;
    }
    
  month--;  
  
  for (i = 0; i < month; i++)
    {
      date += days_in_month[i];
      if (IS_LEAP_YEAR (year) && i == 1) // February of leap year
        date++;
    }
  date += day;
  
  days_in_current_month = days_in_month[month];
  if (IS_LEAP_YEAR (year) && month == 1) // February of leap year
    days_in_current_month++;
  
  
  // check day
  if (day < 1 || day > days_in_current_month)
    {
      fprintf (stderr, "Error: invalid (nonexistent) date\n");
      return -1;
    }
  
  return date;
}
/*!
  \brief Date in format %d.%d.%d reads from text buffer BUF\n
         into date_sim
  \return 0 - in case of error\n
          1 - in case of success.
*/
int
key_read_date (const char *buf, date_sim &date)
{
  int ret_code;
  int day, month, year;
  int days_in_current_month;
  static const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int i;
  
  if (sscanf (buf, "%d.%d.%d", &day, &month, &year) != 3)
    {
      switch (ret_code = key_read_date_ecl (buf, day, month, year))
      {
        case -1:
          {
            fprintf (stderr, "Error: invalid date format: %s\n", buf);
            return -1;
	        }
	      case -4:
	        {
	          fprintf (stderr, "Error: invalid (nonexistent) date: %s\n", buf);
	          return -1;
	        }  
	    }
    }
 
  // Check month and year
  
  if (year < 1900)
    {
      fprintf (stderr, "Error: years earlier 1900 not supported: %s\n", buf);
      return -1;
    }
  
  if (month < 1 || month > 12)
    {
      fprintf (stderr, "Error: invalid (nonexistent) date: %s\n", buf);
      return -1;
    } 
  
  date = 0;
  
  for (i = 1900; i < year; i++)
    {
    if (IS_LEAP_YEAR (i))
      date += 366;
    else
      date += 365;
    }
    
  month--;  
  
  for (i = 0; i < month; i++)
    {
      date += days_in_month[i];
      if (IS_LEAP_YEAR (year) && i == 1) // February of leap year
        date++;
    }
  date += day;
  
  days_in_current_month = days_in_month[month];
  if (IS_LEAP_YEAR (year) && month == 1) // February of leap year
    days_in_current_month++;
  
  
  // check day
  if (day < 1 || day > days_in_current_month)
    {
      fprintf (stderr, "Error: invalid (nonexistent) date: %s\n", buf);
      return -1;
    }
  
  return 0;
}

double 
get_date_day_month_year(date_sim date, int &day, int &month, int &year)
{
  static const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  double next_date;
  
  for (year = 0; date > 0; year++)
    {
      if (IS_LEAP_YEAR (1900 + year))
        next_date = date - 366;
      else
        next_date = date - 365;
      if (next_date <= 0)
        break;
      date = next_date;
    }

  for (month = 0; date > 0; month++)
    {
      next_date = date - days_in_month[month];
      if (IS_LEAP_YEAR (1900 + year) && month == 1)
        next_date--;
      if (next_date <= 0)
        break;
      date = next_date;
    }
  
  day = (int) date;
  month++;
  year += 1900;
    
  return (date - day);
}

void 
d2hms (date_sim date, int &h, int &m, int &s)
{
  double d, n, hh, mm, ss;
  
  d = modf (date, &n);
  d = modf (d * 24.0, &hh);
  d = modf (d * 60.0, &mm);
  d = modf (d * 60.0, &ss);
  h = (int)hh;
  m = (int)mm;
  s = (int)ss;
}

/*!
  \brief Print text from BUF and date from CUR to output 
*/
void
print_date (const date_sim cur, char *buf)
{
  int d, m, y;
  get_date_day_month_year (cur, d, m, y);
  sprintf (buf, "%d.%d.%d", d, m, y);
}

/*!
  \brief print date in ecl format 
*/
void
print_date_ecl (const date_sim cur, char *buf)
{
  int d, m, y;
  get_date_day_month_year(cur, d, m, y);
  
  switch (m)
  {
    case 1:
      sprintf (buf, "%2.2d \'JAN\' %4.4d /", d, y);
      break;
    case 2:
      sprintf (buf, "%2.2d \'FEB\' %4.4d /", d, y);
      break;
    case 3:
      sprintf (buf, "%2.2d \'MAR\' %4.4d /", d, y);
      break;
    case 4:
      sprintf (buf, "%2.2d \'APR\' %4.4d /", d, y);
      break;
    case 5:
      sprintf (buf, "%2.2d \'MAY\' %4.4d /", d, y);
      break;
    case 6:
      sprintf (buf, "%2.2d \'JUN\' %4.4d /", d, y);
      break;
    case 7:
      sprintf (buf, "%2.2d \'JUL\' %4.4d /", d, y);
      break;
    case 8:
      sprintf (buf, "%2.2d \'AUG\' %4.4d /", d, y);
      break;
    case 9:
      sprintf (buf, "%2.2d \'SEP\' %4.4d /", d, y);
      break;
    case 10:
      sprintf (buf, "%2.2d \'OCT\' %4.4d /", d, y);
      break;
    case 11:
      sprintf (buf, "%2.2d \'NOV\' %4.4d /", d, y);
      break;
    case 12:
      sprintf (buf, "%2.2d \'DEC\' %4.4d /", d, y);
      break;
  }
}
