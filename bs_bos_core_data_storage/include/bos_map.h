#ifndef BOS_VAL_TABLE_H
#define BOS_VAL_TABLE_H

namespace blue_sky
  {
  /*!
   \class bos_val_table<T>
   \brief blue_sky map of item_t, accessed with index_t
  */
  template< class index_t, class item_t >
  class bos_val_table : public BS_MAP(index_t, typename dt_val_traits< item_t >::value_type), public objbase
    {
    public:
      //! type of value
      typedef typename dt_val_traits< item_t >::value_type value_type;
      //typedef T value_type;
      //! type of key
      typedef index_t key_type;
      //! type of parent class
      typedef BS_MAP(key_type, value_type) container;
      //! type of map's pair
      typedef typename container::value_type data_pair;
      //! type of reference to second object (value)
      typedef typename container::mapped_type& reference;
      typedef const typename container::mapped_type& const_reference;
      //! type of map's iterator
      typedef typename container::iterator iterator;
      //! type of map's const iterator
      typedef typename container::const_iterator const_iterator;

      using container::erase;
      using container::find;
      using container::begin;
      using container::end;

      //!	Empty destructor.
      virtual ~bos_val_table() {};

      //const_iterator begin () const {
      //	return begin ();
      //}

      //const_iterator end () const {
      //	return end ();
      //}

      /*!
      	\brief Add item method.
      	\param key - key object
      	\param value - value object
      	\return such as std::map
       */
      bool add_item(const key_type& key, const value_type& value)
      {
        return insert(data_pair(key, value)).second;
      }
      /*!
      	\brief Remove item method.
      	\param key - key object
       */
      void rem_item(const key_type& key)
      {
        erase(key);
      }

      /*!
      	\brief Search for item method.
      	\param key - key object
       */
      reference at(const key_type& key)
      {
        return at_internal< iterator, reference >(key);
      }

      const_reference at(const key_type& key) const
        {
          return const_cast< this_t* >(this)->at_internal< const_iterator, const_reference >(key);
        }

    private:
      typedef bos_val_table< index_t, item_t > this_t;

      template< class iterator_t, class ref_t >
      ref_t at_internal(const key_type& key)
      {
        iterator_t p_res(find(key));
        if (p_res != end())
          {
            return p_res->second;
          }
        else
          {
            std::ostringstream ostr;
            ostr << "str_val_table: no element found with key = " << key;
            throw std::out_of_range(ostr.str());
          }
      }

      BLUE_SKY_TYPE_STD_CREATE_T_MEM(bos_val_table)
      BLUE_SKY_TYPE_STD_COPY_T_MEM(bos_val_table)

      BLUE_SKY_TYPE_DECL_T_MEM(bos_val_table, objbase, "bos_val_table",
                               "Array of values of the item_t type indexed by index_t key", "")
    };

  //default ctor implementation
  template< class index_t, class item_t >
  bos_val_table< index_t, item_t >::bos_val_table(bs_type_ctor_param param)
      : objbase(param)
  {}

  //copy ctor implementation
  template< class index_t, class item_t >
  bos_val_table< index_t, item_t >::bos_val_table(const bos_val_table< index_t, item_t >& src) : bs_refcounter (src),
      container(src), objbase(src)
  {}
}

#endif // BOS_VAL_TABLE_H
