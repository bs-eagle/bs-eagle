/*!
   \file bs_array_pool.h
   \brief Contains class for array
*/

#ifndef BS_ARRAY_POOL_H
#define BS_ARRAY_POOL_H

#include "bos_map.h"
#include "throw_exception.h"

#include "shared_array.h"

namespace blue_sky
  {

  enum   //! indexes for dimension parameters
  {
    ARRAY_POOL_NX_A,
    ARRAY_POOL_NX_B,
    ARRAY_POOL_NY_A,
    ARRAY_POOL_NY_B,
    ARRAY_POOL_NZ_A,
    ARRAY_POOL_NZ_B,
    ARRAY_POOL_TOTAL
  };

  /*!
       class bs_array_map< T >
       \brief this class is the map for input data arrays
  */
  template<typename index_type, typename item_type>
  class BS_API_PLUGIN bs_array_map  :  public objbase
  {
    //! Declaration of some objbase methods and members
    BLUE_SKY_TYPE_DECL_T(bs_array_map)

  public:
    //! \brief empty destructor
    ~bs_array_map ();

    typedef index_type                            key_type;       //! type of key
    typedef index_type                            index_t;
    typedef item_type                             item_t;         //! type of value
    typedef bs_array_map <index_type, item_type>  this_t;         //! type of *this
    typedef shared_vector <item_t>                shared_vector_t;

    //! \brief class of array dimension parameters
    struct array_info
    {
      array_info ()
      : def_value(0)
      {
        dimens.assign (0);
      }

      array_info (shared_vector_t &a, const index_t *dim, item_t def_val)
      : array(a)
      , def_value(def_val)
      {
        for (index_t i = 0; i < (index_t)dimens.size (); ++i)
          dimens[i] = dim[i];
      }

      shared_vector_t                               array;
      boost::array <index_t, ARRAY_POOL_TOTAL>      dimens;
      item_t                                        def_value;
    };

    
    typedef bos_val_table <key_type, array_info>  container;            //! type of map class
    typedef typename container::data_pair         value_type;           //! type of value
    typedef smart_ptr<container>                  sp_container_t;       //! type of smart pointer of map
    typedef typename container::iterator          val_table_iterator;   //! type of map's iterator
    typedef typename container::const_iterator    val_table_citerator;  //! type of map's constant iterator
    typedef typename container::reference         val_table_reference;  //! type of map's reference
    typedef typename container::reference         val_table_creference; //! type of map's constant reference
    typedef typename container::size_type         size_type;
    typedef typename container::difference_type   difference_type;


    //! \brief initialization arrays parameters from table of value
    void 
    init(index_t nx, index_t ny, index_t nz);

    //! \brief get array dimension parameters
    void 
    get_dimens (key_type index, index_t &nx1, index_t &ny1, index_t &nz1, index_t &nlen) const
    {
      if (contain (index))
        {
          const boost::array <index_t, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          nx1 = size_p[ARRAY_POOL_NX_A] * ni + size_p[ARRAY_POOL_NX_B];
          ny1 = size_p[ARRAY_POOL_NY_A] * nj + size_p[ARRAY_POOL_NY_B];
          nz1 = size_p[ARRAY_POOL_NZ_A] * nk + size_p[ARRAY_POOL_NZ_B];
          nlen = nx1 * ny1 * nz1;
          return;
        }
      nlen = -1;
    }

    index_t
    get_nx (key_type index) const
    {
      if (contain (index))
        {
          const boost::array <index_t, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          return size_p[ARRAY_POOL_NX_A] * ni + size_p[ARRAY_POOL_NX_B];
        }
      return -1;
    }

    index_t 
    get_ny (key_type index) const
    {
      if (contain (index))
        {
          const boost::array <index_t, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          return size_p[ARRAY_POOL_NY_A] * nj + size_p[ARRAY_POOL_NY_B];
        }
      return -1;
    }

    index_t
    get_nz (key_type index) const
    {
      if (contain (index))
        {
          const boost::array <index_t, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          return size_p[ARRAY_POOL_NZ_A] * nk + size_p[ARRAY_POOL_NZ_B];
        }
      return -1;
    }

    index_t 
    get_nlen(key_type index) const
    {
      if (contain (index))
        {
          const boost::array <index_t, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          index_t nx1 = size_p[ARRAY_POOL_NX_A] * ni + size_p[ARRAY_POOL_NX_B];
          index_t ny1 = size_p[ARRAY_POOL_NY_A] * nj + size_p[ARRAY_POOL_NY_B];
          index_t nz1 = size_p[ARRAY_POOL_NZ_A] * nk + size_p[ARRAY_POOL_NZ_B];
          return nx1 * ny1 * nz1;
        }
      return -1;
    }

    /*!
    	\brief Add item method.
    	\param key - key object
    	\param value - value object
    	\return such as std::map
    */
    bool 
    add_item (const key_type key, shared_vector_t a, const index_t *dimens, item_t def_val)
    {
      bool b = array_map->add_item (key, array_info (a, dimens, def_val));
      if (!b)
        {
          bs_throw_exception ("Can't add array to pool");
        }

      a.assign (def_val);
      return b;
    }

    void 
    create_item (const key_type key, const index_t *dimens, item_t def_val)
    {
      index_t nlen = (dimens[ARRAY_POOL_NX_A] * ni + dimens[ARRAY_POOL_NX_B])
                   * (dimens[ARRAY_POOL_NY_A] * nj + dimens[ARRAY_POOL_NY_B])
                   * (dimens[ARRAY_POOL_NZ_A] * nk + dimens[ARRAY_POOL_NZ_B]);
      item_t *p = 0;
      if (nlen <= 0)
        {
          bs_throw_exception (boost::format (
          "Array length less or equal to zero (pool dimens [%d, %d, %d], array dimens [%d, %d, %d, %d, %d, %d])") \
          % ni % nj % nk % \
          dimens[ARRAY_POOL_NX_A] % dimens[ARRAY_POOL_NX_B] % \
          dimens[ARRAY_POOL_NY_A] % dimens[ARRAY_POOL_NY_B] % \
          dimens[ARRAY_POOL_NZ_A] % dimens[ARRAY_POOL_NZ_B]);
        }
      
      if ((dimens[ARRAY_POOL_NX_A] > 0 && ni == 0) ||
          (dimens[ARRAY_POOL_NY_A] > 0 && nj == 0) ||
          (dimens[ARRAY_POOL_NZ_A] > 0 && nk == 0))
        {
          bs_throw_exception (boost::format (
          "Some used pool dimens are zero (pool dimens [%d, %d, %d], array dimens [%d, %d, %d, %d, %d, %d])") \
          % ni % nj % nk % \
          dimens[ARRAY_POOL_NX_A] % dimens[ARRAY_POOL_NX_B] % \
          dimens[ARRAY_POOL_NY_A] % dimens[ARRAY_POOL_NY_B] % \
          dimens[ARRAY_POOL_NZ_A] % dimens[ARRAY_POOL_NZ_B]);
        }  
        
      //if (contain (key))
      //  {
      //    array_info &info = (*array_map)[key];
      //    if (info.array.size ())
      //      {
      //        p = info.array.data ();
      //      }
      //  }

      //shared_array <item_t> a (allocate (nlen), nlen);
      add_item (key, shared_vector <item_t> (nlen, def_val), dimens, def_val);

      //deallocate (p);
    }

    /*!
    	\brief Remove item method.
    	\param key - value object
    */
    void 
    rem_item (const key_type key)
    {
      array_map->rem_item (key);
    }

    /*!
    	\brief Search for item method.
    	\param key - key object
    */
    val_table_creference 
    at (const key_type key) const
    {
      return array_map->at (key);
    }

    //! \brief begin() to iterate through props arrays in stl-style
    val_table_citerator 
    begin () const
    {
      return array_map->begin ();
    }

    //! \brief end() to iterate through props arrays in stl-style
    val_table_citerator 
    end () const
    {
      return array_map->end ();
    }

    //! \brief size() to get size of arrays map in stl-style
    size_t 
    size () const
    {
      return array_map->size ();
    }

    /*!
    	\brief Overloaded operator [].
    	\return reference to value-object
    */
    val_table_creference 
    operator[] (const key_type &i) const
    {
      return array_map->operator[](i);
    }

    val_table_creference 
    get (const key_type &i) const
    {
      return operator[] (i);
    }

    //! \brief Assignment operator
    this_t &
    operator=(const this_t &src)
    {
      ni = src.ni;
      nj = src.nj;
      nk = src.nk;

      array_map = src.array_map; //give_kernel::Instance().create_object_copy(src.array_map);
      return *this;
    }

    sp_container_t 
    get_array() const
    {
      return array_map;
    }

    bool 
    contain (key_type key) const
    {
      return (array_map->find (key) != array_map->end ());
    }

  private:
    //item_t *
    //allocate (index_t len)
    //{
    //  item_t *p = new item_t [len];
    //  if (!p)
    //    {
    //      bs_throw_exception ("Can't allocate memory");
    //    }

    //  memory_list_.push_back (p);
    //  return p;
    //}

    //void 
    //deallocate (item_t *p)
    //{
    //  free_pointer (p);

    //  for (size_t i = 0, cnt = memory_list_.size (); i < cnt; ++i)
    //    {
    //      if (memory_list_[i] == p)
    //        {
    //          memory_list_.erase (memory_list_.begin () + i);
    //          break;
    //        }
    //    }
    //}

    //void
    //free_pointer (item_t *p)
    //{
    //  delete []p;
    //}

  private:
    
    //std::vector <item_t*>   memory_list_;

    auto_value <index_t, 0> ni, nj, nk;   //! model size
    sp_container_t          array_map;
  };

  namespace tools {

    /** 
     * \brief return array or throw exception if passed array is empty
     * */
    template <typename array_t>
    array_t &
    get_non_empty (array_t &a, const std::string &name="")
    {
      if (a.empty ())
        {
          if (name.length ())
            {
              bs_throw_exception ("Passed array (" + name + ") is empty");
            }
          else
            {
              bs_throw_exception ("Passed array is empty");
            }
        }

      return a;
    }

  } //namespace tools

  bool 
  register_array_map (const plugin_descriptor& pd);
  //void py_export_array_map();
}

#endif

