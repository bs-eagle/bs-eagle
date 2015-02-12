/** 
 * \file export_python_wrapper.h
 * \brief
 * \author Sergey Miryanov
 * \date 08.05.2009
 * */
#ifndef BS_BOS_CORE_BASE_EXPORT_PYTHON_WRAPPER_H_
#define BS_BOS_CORE_BASE_EXPORT_PYTHON_WRAPPER_H_

#include "construct_python_object.h"
#include "python_method_wrapper.h"
#include "python_class_wrapper.h"
#include "bs_object_base.h"

#define PY_EXPORTER(exporter_name_, base_exporter_) \
  template <typename T>                             \
  struct exporter_name_                             \
  {                                                 \
    template <typename class_t>                     \
    static class_t &                                \
    export_class (class_t &class__)                 \
    {                                               \
      using namespace boost::python;                \
                                                    \
      base_exporter_ <T>::export_class (class__)

#define PY_EXPORTER_END                             \
      ;                                             \
      return class__;                               \
    }                                               \
  };

#define GET_BY_VALUE(X) \
  make_getter (X, return_value_policy <return_by_value> ())

namespace blue_sky {
namespace python {

  template <typename T>
  struct default_exporter
  {
    template <typename class_t>
    static class_t &
    export_class (class_t &class__)
    {
      using namespace boost::python;
      class__
        .def ("__init__", make_constructor (&construct_python_object <T>))
        //.def ("__cons__", make_constructor (construct_python_object <T>))
        //.def ("__init__", make_function (init_python_object  <T>))
        ;

      return class__;
    }
  };

  template <typename T>
  struct empty_exporter
  {
    template <typename class_t>
    static class_t &
    export_class (class_t &class__)
    {
      return class__;
    }
  };

  namespace class_type
  {
    template <typename class_t, typename base_t>
    struct abstract_class
    {
      typedef boost::python::class_ <class_t, base_t, boost::noncopyable>  bp_class_t;

      static bp_class_t
      export_class (const std::string &name)
      {
        using namespace boost::python;
        return bp_class_t (name.c_str (), no_init);
      }
      template <typename>
      static void
      register_base_ptr ()
      {
      }
      template <typename, typename>
      static void
      register_class_ptr ()
      {
      }
    };
    template <typename class_t, typename base_t>
    struct concrete_class
    {
      typedef boost::python::class_ <class_t, base_t>  bp_class_t;

      static bp_class_t
      export_class (const std::string &name)
      {
        return bp_class_t (name.c_str ());
      }
      template <typename>
      static void
      register_base_ptr ()
      {
      }
      template <typename, typename>
      static void
      register_class_ptr ()
      {
      }
    };
    template <typename class_t, typename base_t>
    struct blue_sky_class
    {
      typedef boost::python::class_ <class_t, base_t, boost::noncopyable>  bp_class_t;

      static bp_class_t
      export_class (const std::string &name)
      {
        using namespace boost::python;
        return bp_class_t (name.c_str (), no_init);
      }

      template <typename class_t_>
      static void
      register_base_ptr ()
      {
        using namespace boost::python;
        register_ptr_to_python <smart_ptr <class_t_, true> > ();
      }
      template <typename class_t_, typename base_t_>
      static void
      register_class_ptr ()
      {
        using namespace boost::python;
        register_ptr_to_python <smart_ptr <class_t_, true> > ();
        implicitly_convertible <smart_ptr <class_t_, true>, smart_ptr <base_t_, true> > ();
        implicitly_convertible <smart_ptr <objbase, true>, smart_ptr <class_t_, true> > ();
      }
    };
  }

  template <typename class_t, template <typename> class post_exporter = default_exporter, template <typename, typename> class bp_class_type = class_type::concrete_class>
  struct base_exporter_nobs
  {
    typedef bp_class_type <class_t, boost::python::bases <> >   bp_class_type_t;
    typedef typename bp_class_type_t::bp_class_t                bp_class_t;

    static void 
    export_class (const std::string &name)
    {
      using namespace boost::python;

      bp_class_t class__ = bp_class_type_t::export_class (name.c_str ())
        ;

      post_exporter <class_t>::export_class (class__);
      bp_class_type_t::template register_base_ptr <class_t> ();
    }
  };
  template <typename class_t, template <typename> class post_exporter = default_exporter, template <typename, typename> class bp_class_type = class_type::blue_sky_class>
  struct base_exporter
  {
    typedef bp_class_type <class_t, boost::python::bases <objbase> >   bp_class_type_t;
    typedef typename bp_class_type_t::bp_class_t                bp_class_t;

    static void 
    export_class (const std::string &name)
    {
      using namespace boost::python;

      bp_class_t class__ = bp_class_type_t::export_class (name.c_str ())
        ;

      post_exporter <class_t>::export_class (class__);
      bp_class_type_t::template register_base_ptr <class_t> ();
    }
  };
  template <typename class_t, typename base_t, template <typename> class post_exporter = default_exporter, template <typename, typename> class bp_class_type = class_type::blue_sky_class>
  struct class_exporter
  {
    typedef bp_class_type <class_t, boost::python::bases <base_t> >   bp_class_type_t;
    typedef typename bp_class_type_t::bp_class_t                      bp_class_t;

    static void 
    export_class (const std::string &name)
    {
      using namespace boost::python;

      bp_class_t class__ = bp_class_type_t::export_class (name.c_str ())
        ;

      post_exporter <class_t>::export_class (class__);
      bp_class_type_t::template register_class_ptr <class_t, base_t> ();
    }
  };


} // namespace python
} // namespace blue_sky


#endif // #ifndef BS_BOS_CORE_BASE_EXPORT_PYTHON_WRAPPER_H_
