/**
 *       \file  debug_tools.cpp
 *      \brief  Useful tools for debug
 *     \author  Borschuk Oleg
 *       \date  19.06.2006
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */

#ifdef _MPI
#include <mpi.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>

#ifndef UNIX
#include <windows.h>
#include <winbase.h>
#endif //UNIX

#include "debug_tools.h"
#include "debug_macroses.h"

using namespace std;

//! This function summons a MessageBox window with errormessage specified
//! by the parameters of the function and exits the whole program with code -1
/*!
  \brief Error messager
  \param *filename
  \param linenum
  \param *errormessage
*/
int
write_vector_to_file (const char *file_name, double *vector, const int n)
{
  //FILE *fp = 0;
  int i;

  CH_PTR (vector);

  //fp = fopen (file_name, "w");
  ofstream fp(file_name);
  if (!fp)
    //TODO: write error message
    return -1;
  for (i = 0; i < n; ++i)
    //fprintf (fp, "%.20lf\n", vector[i]);
    fp << setprecision(20) << vector[i] << endl;

  //fclose (fp);
  fp.close();

  return 0;
}

int
write_int_vector_to_file (const char *file_name, int* vector, const int n)
{
  //FILE *fp = 0;
  int i;

  CH_PTR (vector);

  //fp = fopen (file_name, "w");
  ofstream fp(file_name);

  if (!fp)
    //TODO: write error message
    return -1;

  for (i = 0; i < n; ++i)
    //fprintf (fp, "%d\n", vector[i]);
    fp << vector[i] << endl;

  //fclose (fp);
  fp.close();

  return 0;
}


