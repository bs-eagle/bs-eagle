/**
 *       \file  py_data_storage_interface.cpp
 *      \brief  Python wrappers data_serializer, data_storage_interface
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */
#ifdef BSPY_EXPORTING_PLUGIN
#include "py_data_storage_interface.h"
#include <boost/python.hpp>

namespace blue_sky
  {
  namespace python
    {

    //////////////////////////////////////////////////////////////////////////
    struct BS_API_PLUGIN data_storage_proxy : public data_storage
      {
        typedef py_data_storage             py_impl_t;
        typedef smart_ptr <py_impl_t>       sp_impl_t;
        typedef data_storage                ds_t;

        BLUE_SKY_TYPE_DECL (data_storage_proxy);

        virtual ~data_storage_proxy () {}

        ds_t &save (const std::string &name, const std::string &value);
        void set_impl (const sp_impl_t &impl);

private:

        sp_impl_t impl_;
      };

    data_storage_proxy::data_storage_proxy (bs_type_ctor_param /*param = NULL */)
    {

    }
    data_storage_proxy::data_storage_proxy (const data_storage_proxy &p)
        : bs_refcounter (p), data_storage (p), impl_ (p.impl_)
    {

    }
    data_storage &
    data_storage_proxy::save (const std::string &name, const std::string &value)
    {
      BS_ASSERT (impl_) (name) (value);

      impl_->save (name, value);
      return *this;
    }
    void
    data_storage_proxy::set_impl (const sp_impl_t &impl)
    {
      impl_ = impl;
    }
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    py_data_storage::py_data_storage (PyObject *self)
        : base_t (data_storage_proxy::bs_type ())
        , self_ (self)
    {
      Py_INCREF (self_);
      get_spx <data_storage_proxy> ()->set_impl (this);
    }

    py_data_storage::~py_data_storage ()
    {
      Py_DECREF (self_);
    }

    void
    py_data_storage::save (const std::string &name, const std::string &value)
    {
      using namespace boost::python;
      call_method <void> (self_, "save", name, value);
    }

    void
    py_data_storage::save_default (const std::string &/*name*/, const std::string &/*value*/)
    {
      BS_ASSERT (false && "PURE CALL");
    }

    void
    py_data_serializer::save (const sp_storage_t &/*storage*/, const sp_obj &/*obj*/) const
    {
      //Sergey Miryanov:
      //  today (23.07.2008) we haven't any way to convert sp_obj
      //  to corresponding python wrapper
      //  in future - if we build table between such types we can
      //  convert python wrappers to objects and vice versa
    }

    //void
    //py_data_storage_interface::register_serializer (py_data_serializer serializer)
    //{
    //  get_lspx (this)->register_serializer (serializer.get_spx <py_data_serializer::wrapped_t> ());
    //}
    void
    py_data_storage_interface::set_storage (const py_data_storage &storage)
    {
      get_spx (this)->set_storage (storage.get_spx <data_storage> ());
    }
    void
    py_data_storage_interface::save (const py_objbase &obj)
    {
      get_spx (this)->save (obj.get_sp ());
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (data_storage_proxy);
    BLUE_SKY_TYPE_STD_COPY (data_storage_proxy);
    BLUE_SKY_TYPE_IMPL (data_storage_proxy, data_storage, "py_data_storage::data_storage_proxy", "py_data_storage::data_storage_proxy", "py_data_storage::data_storage_proxy");

    bool
    data_storage_proxy_register_type (const blue_sky::plugin_descriptor &pd)
    {
      bool res = true;
      res &= BS_KERNEL.register_type (pd, data_storage_proxy::bs_type ());
      BS_ASSERT (res);

      return res;
    }

    //////////////////////////////////////////////////////////////////////////
    void
    py_export_data_storage_interface ()
    {
      using namespace boost::python;

      class_ <py_data_storage, bases <py_objbase>, py_data_storage, boost::noncopyable> ("data_storage"/*, init <PyObject*> ()*/)
      .def ("save", &py_data_storage::save_default)
      ;

      class_ <py_data_storage_interface, bases <py_data_storage_interface::base_t> > ("ds_interface")
      .def ("set_storage", &py_data_storage_interface::set_storage)
      .def ("save", &py_data_storage_interface::save)
      ;
    }

  } // namespace python
} // namespace blue_sky
#endif
