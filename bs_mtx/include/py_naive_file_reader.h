/**
 * \file py_naive_file_reader.h
 * \breif python wrapper for naive_file_reader
 * \author Sergey Miryanov
 * \date 18.06.2008
 * */
#ifndef PY_BS_NAIVE_FILE_READER_H_
#define PY_BS_NAIVE_FILE_READER_H_
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky
  {
  namespace python
    {

    void
    py_export_naive_file_reader ();

  } // namespace python
} // namespace blue_sky


#endif  // #ifdef BSPY_EXPORTING_PLUGIN
#endif  // #ifndef PY_BS_NAIVE_FILE_READER_H_
