/**
 * \file data_storage_interface.h
 * \brief interface for storage different data
 * \author Sergey Miryanov
 * \date 21.07.2008
 * */
#ifndef BS_DATA_STORAGE_INTERFACE_H_
#define BS_DATA_STORAGE_INTERFACE_H_

namespace blue_sky
  {

  class BS_API_PLUGIN data_storage : public objbase
    {
    public:

      virtual ~data_storage () {}
      virtual data_storage &save (const std::string &name, const std::string &value);

      template <typename T>
      data_storage &save (const std::string &name, const T &t)
      {
        return save (name, boost::lexical_cast <std::string> (t));
      }

      BLUE_SKY_TYPE_DECL (data_storage);
    };

  class BS_API_PLUGIN data_serializer : public objbase
    {
    public:

      typedef smart_ptr <data_storage> sp_storage_t;

    public:

      virtual ~data_serializer () {}
      virtual void save (const sp_storage_t &storage, const sp_obj &obj) const;

      const type_descriptor &handled_type () const;

      BLUE_SKY_TYPE_DECL (data_serializer);

    protected:

      type_descriptor handled_type_;
    };

  class BS_API_PLUGIN data_storage_interface : public objbase
    {
    public:

      typedef sp_obj																				sp_storage_t;
      typedef smart_ptr <data_serializer>										sp_serializer_t;
      typedef bos_val_table <std::string, sp_serializer_t>	serializer_list_t;
      typedef smart_ptr <serializer_list_t>									sp_serializer_list_t;

    public:
      virtual ~data_storage_interface () {}

      void save (const sp_obj &obj) const;

      void register_serializer (const sp_serializer_t &serializer);
      void set_storage (const sp_storage_t &storage);

      BLUE_SKY_TYPE_DECL (data_storage_interface);

    private:

      sp_storage_t						storage_;
      sp_serializer_list_t		serializer_list_;
    };

  //////////////////////////////////////////////////////////////////////////
  bool data_storage_register_type (const plugin_descriptor &pd);
  bool data_serializer_register_type (const plugin_descriptor &pd);
  bool data_storage_interface_register_type (const plugin_descriptor &pd);

}	// namespace blue_sky


#endif	// #ifndef BS_DATA_STORAGE_INTERFACE_H_
