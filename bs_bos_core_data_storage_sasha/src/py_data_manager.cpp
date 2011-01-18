#include "bs_bos_core_data_storage_stdafx.h"

#include "py_data_manager.h"
#include "data_manager.h"
#include "py_data_class.h"

using namespace boost::python;

namespace blue_sky
  {
  namespace python
    {

    //template<class strategy_t>
    //py_data_manager<strategy_t>::py_data_manager(const py_event_manager_t &src)
    //: py_objbase(wrapped_t::bs_type ())
    //{
    //  this->template get_lspx <wrapped_t> ()->init(src.get_spx <em_t> ());
    //}

    template<class strategy_t>
    py_data_manager<strategy_t>::py_data_manager(const sp_dm_t &src)
    : py_objbase(src)
    {
    }

    template<class strategy_t>
    py_data_manager<strategy_t>::~py_data_manager() {}

    //template<class strategy_t>
    //void py_data_manager<strategy_t>::read_keyword_file (const std::string &filename, py_event_manager_t &em, py_mesh_rs<strategy_t> &mesh)
    //{
    //  smart_ptr <em_t> sp1 (em.template get_spx <em_t> ());
    //  smart_ptr <mesh_rs <strategy_t> > sp2 (mesh.template get_spx <mesh_rs <strategy_t> >());
    //  this->template get_lspx<wrapped_t> ()->read_keyword_file (filename, sp1, sp2);
    //}

    template<class strategy_t>
    py_idata<strategy_t> py_data_manager<strategy_t>::get_idata () const
      {
        return py_idata<strategy_t> (this->template get_spx<wrapped_t> ()->data);
      }

    template<class strategy_t>
    void py_export_data_manager_t (const char *class_name)
    {

      typedef py_data_manager <strategy_t>	dm_t;

      class_< dm_t, bases<py_objbase> >(class_name, init<const typename dm_t::sp_dm_t &>())
        //.def ("read_keyword_file", &dm_t::read_keyword_file)
        .add_property ("data", &dm_t::get_idata)
        ;
    }

    void py_export_data_manager ()
    {
      py_export_data_manager_t< base_strategy_fif >("data_manager_f");
      py_export_data_manager_t< base_strategy_did >("data_manager_d");
      py_export_data_manager_t< base_strategy_dif >("data_manager_m");
    }

    template py_data_manager<base_strategy_fif>::py_data_manager(const sp_dm_t&);
    template py_data_manager<base_strategy_did>::py_data_manager(const sp_dm_t&);
    template py_data_manager<base_strategy_dif>::py_data_manager(const sp_dm_t&);

  } // namespace python
} // namespace blue_sky
