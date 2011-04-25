#include "stdafx.h"

//#ifdef _MPI
//#include "mpi_csr_comm.h"
//#include "mpi_csr_matrix.h"
//#endif  // #ifdef _MPI

#ifdef _MPI_MY
#include "mpi_mesh_grdecl.h"
#endif

//#include "bs_import_common.h"

#include "jacobian.h"
#include "calc_model.h"

//#include "py_mpi_vector.h"
//#include "py_mpi_csr_matrix.h"

#include "event_manager.h"
#include "py_event_manager.h"

#include "reservoir_simulator.h"
#include "py_reservoir_simulator.h"

#include "py_calc_model.h"

#include "data_storage_interface.h"
#include "py_data_storage_interface.h"

#include "calc_well.h"
#include "py_calc_well.h"
#include "facility_manager.h"

#include "well_controller.h"
#include "well_limit_operation.h"

#include "calc_rho.h"
#include "calc_well_pressure.h"

#include "py_event_base.h"

#include "py_facility_manager.h"
#include "py_reservoir.h"

#include "py_jacobian.h"
#include "default_well.h"
#include "py_well_factory.h"
#include "py_default_wells.h"

#include "fi_params.h"

#include "prepare_fpu.h"

#include "well_results_storage.h"
#include "fip_results_storage.h"

// FIXME:
//#include "csr_ilu_cfl.h"
//#include "py_csr_ilu_cfl_prec.h"

#include <boost/shared_ptr.hpp>




using namespace blue_sky;

#ifdef BSPY_EXPORTING_PLUGIN
using namespace blue_sky::python;
using namespace boost::python;
#endif

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_bos_core", "1.0.0", "Blue Sky BOS Core Plugin", "Blue Sky realization of Black Oil Simulator", "bs_bos_core")

//  struct BS_API_PLUGIN test_y;
//
//  struct BS_API_PLUGIN test_i
//  {
//    virtual void method () = 0;
//
//    virtual ~test_i () {}
//  };
//
//  struct BS_API_PLUGIN test_x : test_i
//  {
//    virtual ~test_x ()
//    {
//      printf ("~test_x\n");
//    }
//    virtual void method ()
//    {
//      printf ("test_x::method\n");
//    }
//
//    virtual void m (test_y *)
//    {
//      printf ("test_x::m\n");
//    }
//  };
//  struct BS_API_PLUGIN test_z : public test_x
//  {
//    virtual ~test_z ()
//    {
//      printf ("~test_z\n");
//    }
//    virtual void method ()
//    {
//      printf ("test_z::method\n");
//    }
//
//    virtual void m (test_y *)
//    {
//      printf ("test_z::m\n");
//    }
//  };
//
//  template <typename T>
//  struct BS_API_PLUGIN test_wrapper : T, wrapper <T>
//  {
//    ~test_wrapper ()
//    {
//      printf ("~test_wrapper\n");
//    }
//
//    test_wrapper (T *t)
//    :t_ (t)
//    {
//
//    }
//
//    void method ()
//    {
//      if (override f = this->get_override ("method"))
//        {
//          f ();
//        }
//      else
//        {
//          t_->method ();
//        }
//    }
//    void method_default ()
//    {
//      printf ("method_default\n");
//    }
//
//    void m (test_y *y)
//    {
//      if (override f = this->get_override ("m"))
//        {
//          f (*y);
//        }
//      else
//        {
//          t_->m (y);
//        }
//    }
//
//    T *t_;
//  };
//
//  template <typename T>
//  static
//  smart_ptr <test_wrapper <T>, false>
//  make_test_x_wrapper ()
//  {
//    printf ("make_test_x_wrapper\n");
//    return new test_wrapper <T> ();
//  }
//
//  template <typename T>
//  static
//  smart_ptr <test_wrapper <T>, false>
//  make_test_x_wrapper_2 (T *i)
//  {
//    printf ("make_test_x_wrapper_2\n");
//    return new test_wrapper <T> (i);
//  }
//
//  struct BS_API_PLUGIN test_y
//  {
//    virtual ~test_y ()
//    {
//      printf ("~test_y\n");
//    }
//    virtual void method (test_x *x)
//    {
//      printf ("test_y_x::method\n");
//      x->method ();
//    }
//    virtual void mz (test_z *x)
//    {
//      printf ("test_y_z::method\n");
//      x->method ();
//    }
//    virtual void mi (test_i *x)
//    {
//      printf ("test_y_i::method\n");
//      x->method ();
//    }
//    virtual void m2 (const smart_ptr <test_wrapper <test_i>, false> &x)
//    {
//      printf ("test_y::m2\n");
//      x->method ();
//    }
//
//    virtual void m3 (smart_ptr <test_wrapper <test_x>, false> &x)
//    {
//      printf ("test_y::m3\n");
//      x->method ();
//    }
//  };
//
//
//  template <typename T> inline T * get_pointer(smart_ptr <T, false> const & p)
//  {
//    return const_cast <T *> (p.get());
//  }
//
//  using namespace boost::python;
//
//  template <typename T>
//  object init__ (object py_obj)
//  {
//    printf ("init__\n");
//    object return_value = call_method<object>(py_obj.ptr(), "__cons__");
//    initialize_wrapper (py_obj.ptr (), extract <test_wrapper <T> *> (py_obj));
//    return return_value;
//  }
//
//  template <typename T>
//  object init__2 (object py_obj, T *i)
//  {
//    printf ("init__\n");
//    object return_value = call_method<object>(py_obj.ptr(), "__cons__", i);
//    initialize_wrapper (py_obj.ptr (), extract <test_wrapper <T> *> (py_obj));
//    return return_value;
//  }
//
//  void
//  py_export_test_x ()
//  {
//    class_ <test_i, boost::noncopyable> ("test_i", no_init)
//      ;
//
//    class_ <test_x/*, test_wrapper <test_x>*/, bases <test_i>/*, boost::noncopyable*/> ("test_x"/*, boost::python::no_init*/)
//      .def ("method", &test_x::method)//, &test_wrapper<test_x>::method_default)
//      .def ("m", &test_x::m)//, &test_wrapper<test_x>::method_default)
//      //.def ("__cons__", make_constructor (make_test_x_wrapper <test_x>))
//      //.def ("__init__", init__ <test_x>)
//      ;
//    class_ <test_z/*, test_wrapper <test_z>*/, bases </*test_wrapper <*/test_x/*>*/ >/*, boost::noncopyable*/> ("test_z"/*, boost::python::no_init*/)
//      //.def ("__cons__", make_constructor (make_test_x_wrapper <test_z>))
//      //.def ("__init__", init__ <test_z>)
////      .def ("m", &test_z::method)//, &test_wrapper <test_z>::method_default)
//      //.def ("method", &test_z::method)//, &test_wrapper <test_z>::method_default)
//      ;
//
//    class_ <test_wrapper <test_x>, bases <test_i>, boost::noncopyable> ("test_x_wrapper", boost::python::no_init)
//      .def ("__cons__", make_constructor (make_test_x_wrapper_2 <test_x>))
//      .def ("__init__", init__2 <test_x>)
//      .def ("method", &test_wrapper <test_x>::method)//, &test_wrapper <test_z>::method_default)
//      .def ("m", &test_wrapper <test_x>::m)//, &test_wrapper <test_z>::method_default)
//      ;
//
//    class_ <test_wrapper <test_z>, bases <test_x>, boost::noncopyable> ("test_z_wrapper", boost::python::no_init)
//      .def ("__cons__", make_constructor (make_test_x_wrapper_2 <test_z>))
//      .def ("__init__", init__2 <test_z>)
//      .def ("method", &test_wrapper <test_z>::method)//, &test_wrapper <test_z>::method_default)
//      .def ("m", &test_wrapper <test_z>::m)//, &test_wrapper <test_z>::method_default)
//      ;
//
//    class_ <test_y, boost::noncopyable> ("test_y")
//      .def ("method", &test_y::method)
//      //.def ("method", &test_y::mz)
//      .def ("mz", &test_y::mz)
//      .def ("mi", &test_y::mi)
//      .def ("m2", &test_y::m2)
//      .def ("m3", &test_y::m3)
//      ;
//
//    typedef smart_ptr <test_x, false> test_x_t;
//    typedef smart_ptr <test_wrapper <test_x>, false> test_x_wrapper_t;
//    register_ptr_to_python<test_x_t>();
//    register_ptr_to_python<test_x_wrapper_t>();
//
//    typedef smart_ptr <test_z, false> test_z_t;
//    typedef smart_ptr <test_wrapper <test_z>, false> test_z_wrapper_t;
//    register_ptr_to_python<test_z_t>();
//    register_ptr_to_python<test_z_wrapper_t>();
//
//    register_ptr_to_python<test_y *>();
//    register_ptr_to_python<test_i *>();
//    register_ptr_to_python<smart_ptr <test_i, false> >();
//  }

  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    //////////////////////////////Events/////////////////////////////////////////////
    res &= event_base_register_types (pd);
    res &= event_manager_register_types (pd);

    res &= blue_sky::well_events_register_type (pd);
    BS_ASSERT (res);

    //////////////////////////////Linear Solver///////////////////////////////
// FIXME:
//#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
//    res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, csr_ilu_cfl_prec<base_strategy_fi>::bs_type());
//    BS_ASSERT (res);
//    res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, csr_ilu_cfl_prec<base_strategy_di>::bs_type());
//    BS_ASSERT (res);
//    res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, csr_ilu_cfl_prec<base_strategy_mixi>::bs_type());
//    BS_ASSERT (res);
//#endif


//#ifdef _MPI
//    res &= blue_sky::mpi_csr_matrix_register_type (pd);
//    BS_ASSERT (res);
//    res &= blue_sky::mpi_csr_comm_register_type (pd);
//    BS_ASSERT (res);
//#endif

    res &= BS_KERNEL.register_type (pd, well_results_storage::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, fip_results_storage::bs_type ()); BS_ASSERT (res);

    res &= calc_rho_register_types (pd);            BS_ASSERT (res);
    res &= calc_well_pressure_register_types (pd);  BS_ASSERT (res);
    res &= calc_well_register_types (pd);           BS_ASSERT (res);
    res &= wells::default_well_register_types (pd); BS_ASSERT (res);

    res &= facility_manager_register_type (pd);
    BS_ASSERT (res);
    res &= well_factory_register_type (pd);
    BS_ASSERT (res);
    res &= wells::well_controller_factory_register_type (pd);
    BS_ASSERT (res);
    res &= wells::well_limit_operation_factory_register_type (pd);
    BS_ASSERT (res);

	#ifdef BSPY_EXPORTING_PLUGIN
    res &= data_storage_proxy_register_type (pd);
    BS_ASSERT (res);
	#endif

    res &= data_storage_register_type (pd);
    BS_ASSERT (res);
    res &= data_serializer_register_type (pd);
    BS_ASSERT (res);
    res &= data_storage_interface_register_type (pd);
    BS_ASSERT (res);

    res &= reservoir_simulator_register_types (pd);
    BS_ASSERT (res);

    return res;
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {

    //bool res = true;
    return register_types (*bs_init.pd_);

  }
}//bs
//
//namespace blue_sky {
//  template <typename bcsr_d, typename bcsr_f>
//  void
//  convert_double_to_float (bcsr_d *bcsr_mx_d, bcsr_f *bcsr_mx_f)
//  {
//    //typedef typename mx_t::i_type_t i_type;
//    //typedef typename mx_t::fp_type_t fp_type;
//
//    typedef typename bcsr_d::index_t index_t;
//    typedef typename bcsr_d::item_t item_t;
//    typedef typename bcsr_d::item_array_t item_array_t;
//    typedef typename bcsr_d::index_array_t index_array_t;
//
//    //typedef typename bcsr_d::index_t index_t;
//    typedef typename bcsr_f::item_t rhs_item_t;
//
//    const index_t *rows_d = &bcsr_mx_d->get_rows_ptr ()[0];
//    const index_t *cols_d = &bcsr_mx_d->get_cols_ind ()[0];
//    const item_t *val_d   = &bcsr_mx_d->get_values ()[0];
//
//    index_t *rows_f = &bcsr_mx_f->get_rows_ptr ()[0];
//    index_t *cols_f = &bcsr_mx_f->get_cols_ind ()[0];
//    rhs_item_t *val_f   = &bcsr_mx_f->get_values ()[0];
//
//    typedef typename bcsr_d::i_type_t i_type;
//
//    bcsr_mx_f->alloc_rows_ptr (bcsr_mx_d->n_rows, 0);
//    bcsr_mx_f->alloc_cols_ind_and_values (bcsr_mx_d->get_n_non_zeros ());
//
//    bcsr_mx_f->get_rows_ptr ().assign (bcsr_mx_d->get_rows_ptr ().begin (), bcsr_mx_d->get_rows_ptr ().end ());
//    bcsr_mx_f->get_diag_ind ().assign (bcsr_mx_d->get_diag_ind ().begin (), bcsr_mx_d->get_diag_ind ().end ());
//    bcsr_mx_f->get_cols_ind ().assign (bcsr_mx_d->get_cols_ind ().begin (), bcsr_mx_d->get_cols_ind ().end ());
//    bcsr_mx_f->get_values ().assign (bcsr_mx_d->get_values ().begin (), bcsr_mx_d->get_values ().end ());
//
//    bcsr_mx_f->n_rows = bcsr_mx_d->n_rows;
//    bcsr_mx_f->n_cols = bcsr_mx_d->n_cols;
//    bcsr_mx_f->n_block_size = bcsr_mx_d->n_block_size;
//
//    index_array_t tmp_rows_ptr, tmp_diag_ind, tmp_cols_ind;
//    item_array_t tmp_values;
//    std::swap (bcsr_mx_d->get_rows_ptr (), tmp_rows_ptr);
//    std::swap (bcsr_mx_d->get_diag_ind (), tmp_diag_ind);
//    std::swap (bcsr_mx_d->get_cols_ind (), tmp_cols_ind);
//    std::swap (bcsr_mx_d->get_values (), tmp_values);
//  }
//
//  template <typename py_bcsr_d, typename py_bcsr_f>
//  void
//  py_convert_double_to_float (py_bcsr_d &bcsr_d, py_bcsr_f &bcsr_f)
//  {
//    convert_double_to_float (&(*bcsr_d.template get_lspx <typename py_bcsr_d::matrix_t> ()),
//      &(*bcsr_f.template get_lspx <typename py_bcsr_f::matrix_t> ()));
//  }
//}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  py_export_fi_params ();

//#ifdef _MPI
//  python::py_export_mpi_vector ();
//  python::py_export_mpi_csr_matrix ();
//#endif


  py_export_events ();
  py_export_event_manager();

  python::py_export_reservoir_simulator ();
  python::py_export_calc_model ();
  // we should export enum only for one instance of reservoir simulator
  //reservoir_simulator<base_strategy_di>::py_export_signals_enum (); //!TODO:

  python::py_export_calc_well ();

  python::py_export_data_storage_interface ();

  python::py_export_jacobian ();

  python::py_export_facility_manager ();
  python::py_export_reservoir ();

  python::py_export_well_factories ();

// FIXME:
//#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
//  python::py_export_csr_ilu_cfl_prec ();
//#endif

  python::py_export_default_wells ();

  def ("enable_fpu_exceptions", blue_sky::tools::prepare_fpu::enable_exceptions);

  //def ("convert_double_to_float", py_convert_double_to_float <python::py_bcsr_matrix <base_strategy_di::item_array_t, base_strategy_di::index_array_t>, python::py_bcsr_matrix <base_strategy_fi::item_array_t, base_strategy_fi::index_array_t> >);

  //py_export_test_x ();
}
#endif //BSPY_EXPORT_PLUGIN
