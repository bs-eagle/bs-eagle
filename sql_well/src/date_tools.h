#ifndef DATE_TOOLS_H
#define DATE_TOOLS_H

/*!
  \file date_tools.h
  \brief Declarations of date handling tools
*/

#include <time.h>

typedef double date_sim;


int key_read_date (const char *buf, date_sim &date);

int key_read_time (const char *buf, date_sim &tm);

//the function for reading date as (1 jan 1900)
int key_read_date_as_text (const char *buf, int &day, int &month, int &year);


double get_date_day_month_year (date_sim date, int &day, int &month, int &year);


void print_date (const date_sim cur, char *buf);

void print_date_ecl (const date_sim cur, char *buf);


#endif // DATE_TOOLS_H
