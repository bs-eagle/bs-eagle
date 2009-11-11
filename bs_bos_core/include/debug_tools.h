/**
 *       \file  debug_tools.h
 *      \brief  Useful tools for debug
 *     \author  Borschuk Oleg
 *       \date  19.06.2006
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef __DEBUG_TOOLS_H_
#define __DEBUG_TOOLS_H_

#ifdef _MPI
#include <mpi.h>
#endif

#include <stdio.h>

// Generate a messagebox with error message and exit a program
int error_messager (const char *filenum, const int line, const char *errormessage);

/** вспомогательная структура для автоматического закрытия файла */
struct AutoFile
  {
    /** файл для закрытия */
    FILE *fp_;

    AutoFile(FILE *fp = 0) : fp_(fp) {}

    ~AutoFile()
    {
      if (fp_)
        fclose(fp_);
    }

    FILE *operator()() const
      {
        return fp_;
      }
  };

bool file_exists(const char *file_name);

/*!
  \param offset number of items from beginning of file
  \return 0 if successful
 */
int write_vector_to_file(FILE *fp, double *vector, size_t n, int offset);

int write_vector_to_file(const char *file_name, double *vector, size_t n, int offset);

/*!
  \param offset number of items from beginning of file
  \return 0 if successful
 */
int read_vector_from_file(FILE *fp, double *vector, size_t n, int offset);

int read_vector_from_file(const char *file_name, double *vector, size_t n, int offset);

// print vector to file (one string is one value)
// DANGER! overflow is not checked!

int write_vector_to_file (const char *file_name, double *vector, const int n);

// print int vector to file (one string is one value)
// DANGER! overflow is not checked!
int write_int_vector_to_file (const char *file_name, int* vector, const int n);
#ifdef _MPI
int mpi_write_int_vector_to_file (const char *file_name, int* vector, const int n, MPI_Comm comm);
int mpi_write_vector_to_file (const char *file_name, double* vector, const int n, MPI_Comm comm);
#endif
// read vector from file
// DANGER! overflow is not checked!
int read_vector_from_file (const char *file_name, double* vector, const int n);

int write_int_vector_to_file (const char *file_name, int* vector, const int n);
#endif //__DEBUG_TOOLS_H_

