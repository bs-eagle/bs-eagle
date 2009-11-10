/**
 * \file py_save_seq_vector.h
 * \breif python wrapper for save_seq_vector
 * \author Sergey Miryanov
 * \date 20.06.2008
 * */
#ifndef PY_BS_SAVE_SEQ_VECTOR_H_
#define PY_BS_SAVE_SEQ_VECTOR_H_
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky
  {
  namespace python
    {

    void
    py_export_save_seq_vector ();

  } // namespace python
} // namespace blue_sky


#endif  // #ifdef BSPY_EXPORTING_PLUGIN
#endif  // #ifndef PY_BS_SAVE_SEQ_VECTOR_H_
