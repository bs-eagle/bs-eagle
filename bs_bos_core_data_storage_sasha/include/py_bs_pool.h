//#ifndef PY_BS_POOL_H
//#define PY_BS_POOL_H
///**
// * \file   py_bs_pool.h
// * \brief  Python wrapper for template class bs_array_map<value_type>
// * \author Nikonov Maxim
// * \date 2008-05-30
// */
//
//#include "pool.h"
//
//namespace blue_sky
//  {
//  namespace python
//    {
//    template <typename index_t, typename item_array_t>
//    class py_array_map : public py_objbase
//      {
//      public:
//
//        typedef py_objbase                                  based_t;                    //! type of python base class wrapped objbase
//        typedef bs_array_map <index_t, item_array_t>        wrapped_t;                  //! type of wrapped object
//        typedef typename wrapped_t::sp_array_t              sp_array_t;
//        typedef typename wrapped_t::value_type              value_type;                 //! type of value
//        typedef typename wrapped_t::key_type                key_type;                   //! type of key value
//        typedef typename wrapped_t::val_table_creference    val_table_creference;
//        typedef typename wrapped_t::val_table_citerator     val_table_citerator;
//        //typedef typename wrapped_t::second_type             second_type;
//        typedef typename wrapped_t::size_type               size_type;
//        typedef typename wrapped_t::difference_type         difference_type;
//
//        //using wrapped_t::erase;
//
//        py_array_map()
//            : based_t(wrapped_t::bs_type())
//        {}
//
//        void init(int nx, int ny, int nz)
//        {
//          this->get_spx (this)->init(nx,ny,nz);
//        }
//
//        int get_nx(const int index) const
//          {
//            return this->get_spx(this)->get_nx(index);
//          }
//
//        int get_ny(const int index) const
//          {
//            return this->get_spx(this)->get_ny(index);
//          }
//
//        int get_nz(const int index) const
//          {
//            return this->get_spx(this)->get_nz(index);
//          }
//
//        int get_nlen(const int index) const
//          {
//            return this->get_spx(this)->get_nlen(index);
//          }
//
//        bool add_item (const key_type key, const sp_array_t& value, const int *dimens, typename wrapped_t::item_t def_val)
//        {
//          return this->get_spx(this)->add_item(key,value,dimens,def_val);
//        }
//
//        void rem_item (const key_type key)
//        {
//          this->get_spx(this)->rem_item(key);
//        }
//        val_table_creference at (const key_type key) const;
//        val_table_citerator begin () const;
//        val_table_citerator end () const;
//        size_t size () const
//          {
//            return this->get_spx(this)->size();
//          }
//
//        val_table_creference operator[] (const key_type &i);
//
//      };
//
//    template <typename T>
//    void py_export_array_map(const char *class_name)
//    {
//      boost::python::class_<T, boost::python::bases<py_objbase> >(class_name)
//      .def("init",      &T::init)
//      .def("get_nx",    &T::get_nx)
//      .def("get_ny",    &T::get_ny)
//      .def("get_nz",    &T::get_nz)
//      .def("get_nlen",  &T::get_nlen)
//      .def("add_item",  &T::add_item)
//      .def("rem_item",  &T::rem_item)
//      .def("size",      &T::size)
//      ;
//    }
//
//    void py_export_array_maps();
//  } // namespace python
//} // namespace blue_sky
//
//#endif // PY_BS_POOL_H
