/**
 *       \file  py_data_storage_interface.h
 *      \brief  Python wrappers data_serializer, data_storage_interface
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */
#ifndef BS_PY_DATA_STORAGE_INTERFACE_H_
#define BS_PY_DATA_STORAGE_INTERFACE_H_

#include "data_storage_interface.h"
#include "py_bs_object_base.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {


    class BS_API_PLUGIN py_data_storage : public py_objbase
      {
      public:

        typedef py_objbase    base_t;

      public:

        py_data_storage (PyObject *self);
        virtual ~py_data_storage ();

        void save (const std::string &name, const std::string &value);
        void save_default (const std::string &name, const std::string &value);

      private:

        PyObject *self_;
      };

    class BS_API_PLUGIN py_data_serializer : public py_objbase, public data_serializer
      {
      public:

        typedef data_serializer   wrapped_t;
        typedef py_objbase        base_t;

      public:

        py_data_serializer (PyObject *self)
            : base_t (this)
            , self_ (self)
        {

        }

        void save (const sp_storage_t &storage, const sp_obj &obj) const;
        void save_impl (py_data_storage storage, py_objbase obj);

      private:

        PyObject *self_;
      };


    class BS_API_PLUGIN py_data_storage_interface : public py_objbase
      {
      public:

        typedef data_storage_interface  wrapped_t;
        typedef py_objbase              base_t;

      public:

        py_data_storage_interface ()
            : base_t (wrapped_t::bs_type ())
        {

        }


        void save (const py_objbase &obj);

        // see py_data_serializer::save for explanation
        //void register_serializer (py_data_serializer *serializer);

        void set_storage (const py_data_storage &storage);
      };

    bool
    data_storage_proxy_register_type (const blue_sky::plugin_descriptor &pd);

    void
    py_export_data_storage_interface ();

  } // namespace python
} // namespace blue_sky


#endif
#endif  // #ifndef BS_PY_DATA_STORAGE_INTERFACE_H_
