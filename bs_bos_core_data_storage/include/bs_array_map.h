/*!
   \file bs_array_bs_array_map.h
   \brief Contains class for array
*/

#ifndef BS_ARRAY_MAP_H
#define BS_ARRAY_MAP_H

#include "bos_map.h"
#include "throw_exception.h"

#include "auto_value.h"
#include <boost/array.hpp>

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

    typedef std::string                           key_type;       //! type of key
    typedef index_type                            index_t;
    typedef item_type                             item_t;         //! type of value
    typedef bs_array_map <index_type, item_type>  this_t;         //! type of *this
    typedef smart_ptr<bs_array <item_t>, true>    sp_bs_array_t;

    //! \brief class of array dimension parameters
    struct array_info
    {
      array_info ()
      : def_value(0)
      {
        dimens.assign (0);
      }

      array_info (sp_bs_array_t &a, const int *dim, item_t def_val)
      : array(a)
      , def_value(def_val)
      {
        for (index_t i = 0; i < (index_t)dimens.size (); ++i)
          dimens[i] = dim[i];
      }

      sp_bs_array_t                                 array;
      boost::array <int, ARRAY_POOL_TOTAL>          dimens;
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
          const boost::array <int, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
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
          const boost::array <int, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          return size_p[ARRAY_POOL_NX_A] * ni + size_p[ARRAY_POOL_NX_B];
        }
      return -1;
    }

    index_t 
    get_ny (key_type index) const
    {
      if (contain (index))
        {
          const boost::array <int, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          return size_p[ARRAY_POOL_NY_A] * nj + size_p[ARRAY_POOL_NY_B];
        }
      return -1;
    }

    index_t
    get_nz (key_type index) const
    {
      if (contain (index))
        {
          const boost::array <int, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          return size_p[ARRAY_POOL_NZ_A] * nk + size_p[ARRAY_POOL_NZ_B];
        }
      return -1;
    }

    index_t 
    get_nlen(key_type index) const
    {
      if (contain (index))
        {
          const boost::array <int, ARRAY_POOL_TOTAL> &size_p = (*array_map)[index].dimens;
          index_t nx1 = size_p[ARRAY_POOL_NX_A] * ni + size_p[ARRAY_POOL_NX_B];
          index_t ny1 = size_p[ARRAY_POOL_NY_A] * nj + size_p[ARRAY_POOL_NY_B];
          index_t nz1 = size_p[ARRAY_POOL_NZ_A] * nk + size_p[ARRAY_POOL_NZ_B];
          return nx1 * ny1 * nz1;
        }
      return -1;
    }

   
    boost::python::list 
    py_list_items ()
    {
      val_table_iterator it;
      boost::python::list items;
      for (it = array_map->begin(); it != array_map->end(); it++)
        items.append(it->first);
      return items;
    }
    
    
     /*!
    	\brief Add item method.
    	\param key - key object
    	\param value - value object
    	\return such as std::map
    */
    bool 
    add_item (const key_type key, sp_bs_array_t &a, const int *dimens, item_t def_val)
    {
      if (!contain (key))
        {
          if (!array_map->add_item (key, array_info (a, dimens, def_val)))
            {
              bs_throw_exception ("Can't add array to pool");
              return false;
            }

          a->assign (def_val);
          return true;
        }
      else
        {
          bs_throw_exception ("Array " + key + " already exists");
          return false;
        }  
    }

    void 
    py_create_item (const key_type key, boost::python::list dimens, item_t def_val)
      {
        int new_dimens[ARRAY_POOL_TOTAL];
        for (int i = 0; i < ARRAY_POOL_TOTAL; i++) 
          {
            new_dimens[i] = boost::python::extract<int>(dimens[i]);
          }
        create_item (key, &new_dimens[0],def_val);
    }

    sp_bs_array_t 
    create_item (const key_type key, const int *dimens, item_t def_val)
    {
      index_t nlen = (dimens[ARRAY_POOL_NX_A] * ni + dimens[ARRAY_POOL_NX_B])
                   * (dimens[ARRAY_POOL_NY_A] * nj + dimens[ARRAY_POOL_NY_B])
                   * (dimens[ARRAY_POOL_NZ_A] * nk + dimens[ARRAY_POOL_NZ_B]);
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
      sp_bs_array_t a = give_kernel::Instance().create_object(bs_array<item_t>::bs_type());
      a->resize (nlen);
      a->assign (def_val);
       
      add_item (key, a, dimens, def_val);
      return a;
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

  private:
    
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
      if (a->empty ())
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

#endif //BS_ARRAY_MAP_H
