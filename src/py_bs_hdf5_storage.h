/**
 * \file py_bs_hdf5_storage.h
 * \brief Python wrappers for bs_hdf5_storage
 * \author
 * \date
 */

#ifndef BS_PY_BS_HDF5_STORAGE_H_
#define BS_PY_BS_HDF5_STORAGE_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include "py_bs_object_base.h"
#include "hdf5_functions.h"
#include "aligned_allocator.h"
#include "seq_vector.h"
#include <vector>

namespace blue_sky
  {
  namespace python
    {

    class py_bs_hdf5_storage : public py_objbase
      {
      public:
        typedef bs_hdf5_storage wrapped_t;
        //! constructor
        py_bs_hdf5_storage ()
            : py_objbase (BS_KERNEL.create_object (bs_hdf5_storage::bs_type ()))
        {
        }

        int open (const char *filename)
        {
          return get_spx (this)->open (filename);
        }
        int close ()
        {
          return get_spx (this)->close ();
        }

        template<typename T>
        int write (seq_vector<T> &x, std::string dataset_name)
        {
          return write_data_to_hdf5 (get_spx (this)->get_file_id(), dataset_name, &x[0], (int)x.size ());
        }

      };

    //! export wrapper into python
    void py_export_bs_hdf5_storage ();

  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef BS_PY_BS_HDF5_STORAGE_H_

