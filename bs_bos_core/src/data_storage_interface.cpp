/**
 *       \file  data_storage_interface.cpp
 *      \brief  Implementation of data_storage_interface
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "data_storage_interface.h"

#ifdef _DEBUG
#include "well_serializer.h"
#endif


namespace blue_sky
  {

  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bos_val_table <std::string, smart_ptr <data_serializer> >), 1, (objbase), "bos_map <str, data_serializer>", "bos_map <str, data_serializer>", "bos_map <str, data_serializer>", false);

  data_storage::data_storage (bs_type_ctor_param /*param = NULL */)
  {

  }
  data_storage::data_storage (const data_storage & /*ds*/)
        : bs_refcounter (), objbase ()
  {

  }

  data_storage &
  data_storage::save (const std::string &name, const std::string &value)
  {
    BS_ASSERT (false && "BASE METHOD CALL") (name) (value);
    return *this;
  }


  data_serializer::data_serializer (bs_type_ctor_param /*param = NULL */)
  {

  }
  data_serializer::data_serializer (const data_serializer& ds)
        : bs_refcounter (), objbase ()
  {
    handled_type_ = ds.handled_type ();
  }

  void
  data_serializer::save (const sp_storage_t &storage, const sp_obj &obj) const
    {
      BS_ASSERT (storage);
      BS_ASSERT (obj);

      const sp_storage_t &locked_storage (storage);
      data_storage &st = *locked_storage;

      BS_ASSERT (false && "BASE METHOD CALL") (storage->bs_resolve_type ().stype_) (obj->bs_resolve_type ().stype_);
      st.save ("base method call", "---");
    }

  const type_descriptor &
  data_serializer::handled_type () const
    {
      return handled_type_;
    }


  data_storage_interface::data_storage_interface (bs_type_ctor_param /*param = NULL */)
  {
    serializer_list_ = BS_KERNEL.create_object (serializer_list_t::bs_type ());
    BS_ASSERT (serializer_list_);

#ifdef _DEBUG
    register_serializer (BS_KERNEL.create_object (wells::well_serializer::bs_type ()));
#endif
  }
  data_storage_interface::data_storage_interface (const data_storage_interface& dsi)
        : bs_refcounter (), objbase ()
  {
    storage_ = dsi.storage_;
    serializer_list_ = dsi.serializer_list_;

    BS_ASSERT (serializer_list_);
  }

  void
  data_storage_interface::save (const sp_obj &obj) const
    {
      BS_ASSERT (storage_);
      if (!storage_)
        return ;

      BS_ASSERT (obj);
      const type_descriptor &td = obj->bs_resolve_type ();

      serializer_list_t::const_iterator it = serializer_list_->find (td.stype_);
      BS_ASSERT (it != serializer_list_->end ()) (td.stype_);

      if (it != serializer_list_->end ())
        it->second->save (storage_, obj);
    }

  void
  data_storage_interface::register_serializer (const sp_serializer_t &serializer)
  {
    BS_ASSERT (serializer);
    BS_ASSERT (serializer->bs_resolve_type ().stype_.length ());

    const type_descriptor &td = serializer->handled_type ();
    BS_ASSERT (td.stype_.length ()) (serializer->bs_resolve_type ().stype_);

    const sp_serializer_list_t &locked_list (serializer_list_);
    serializer_list_t::iterator it = locked_list->find (td.stype_);
    BS_ASSERT (it == locked_list->end ());

    if (it == locked_list->end ())
      locked_list->insert (make_pair (td.stype_, serializer));
    else
      it->second = serializer;
  }

  void
  data_storage_interface::set_storage (const sp_storage_t &storage)
  {
    storage_ = storage;

    BS_ASSERT (storage_);
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (data_storage);
  BLUE_SKY_TYPE_STD_COPY (data_storage);
  BLUE_SKY_TYPE_IMPL (data_storage, objbase, "data_storage", "data_storage", "data_storage");

  BLUE_SKY_TYPE_STD_CREATE (data_serializer);
  BLUE_SKY_TYPE_STD_COPY (data_serializer);
  BLUE_SKY_TYPE_IMPL (data_serializer, objbase, "data_serializer", "data_serializer", "data_serializer");

  BLUE_SKY_TYPE_STD_CREATE (data_storage_interface);
  BLUE_SKY_TYPE_STD_COPY (data_storage_interface);
  BLUE_SKY_TYPE_IMPL (data_storage_interface, objbase, "data_storage_interface", "data_storage_interface", "data_storage_interface");

  //////////////////////////////////////////////////////////////////////////
  bool
  data_storage_register_type (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, data_storage::bs_type ());
    BS_ASSERT (res);

    return res;
  }
  bool
  data_serializer_register_type (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, data_serializer::bs_type ());
    BS_ASSERT (res);

    return res;
  }
  bool
  data_storage_interface_register_type (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, data_storage_interface::bs_type ());
    BS_ASSERT (res);

    return res;
  }
} // namespace blue_sky
