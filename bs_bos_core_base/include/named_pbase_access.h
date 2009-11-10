/*! \file named_pbase_access.h
		\brief This file contains named_pbase class - for
					 named access to property_base tables
		\author Nikonov Max
*/

#ifndef NAMED_PBASE_ACCESS_H
#define NAMED_PBASE_ACCESS_H

#include "property_base.h"
#include "throw_exception.h"
#include "construct_python_object.h"

/* How to use this:
	 ~~~~~~~~~~~~~~~~
	If you want to inherite from property_base and have named access
	to property table elements, just write:

	SIMPLE C++ CLASS:

		class A : public named_pbase {
		public:
			//blue-sky class declaration
			BLUE_SKY_TYPE_DECL(A);
		public:
			// enumeration definitions
			PROP_BASE_IDX_DECL_BEGIN(A,property_base) //or named_pbase in second argument place
				//your indexes divided with commas
				A1,
				A2,
				A_TOTAL,
			PROP_BASE_IDX_DECL_END

    virtual ~A (); // virtual dtor

		PBASE_ACCESS_MS(A) // method for access property table
		// It defines the idx_type type and uses it to access
		// property table elements

		public:
			void set_default_values (); // Set default values (may be unnecessary)

			virtual const std::string &get_params_name (idx_type idx); // returns name by index
		};

	As implementations write:

		void A::set_default_values () {
			this->params_names.resize(A_TOTAL);
			// result will be for idx_type A1 sets name "A1" and description "A1 descr"
			EREG(A,A1,"A1 descr");
			EREG(A,A2,"A2 special descr");

			resize(AMG_TOTAL); // to resize names table and table of properties

			set_param(A1,(float)0.25); // this instruction sets up parameter A1 in property_table
			set_param(A2,"Hello, World!");
	  }

	const std::string & A::get_params_name (idx_type idx) {
		return params_names[idx].first;
	}

	// And special blue-sky objbase class implementations

	A::A (bs_type_ctor_param param)
		: bs_refcounter(), named_pbase()
	{
		set_default_values ();
	}

	A::A (const A& prop)
		: bs_refcounter(), named_pbase(prop)
	{
		if (&prop != this)
			*this = prop;

		set_default_values ();
	}

	A::~A () {
		set_default_values ();
	}

	BLUE_SKY_TYPE_STD_CREATE(A)
  BLUE_SKY_TYPE_STD_COPY(A)
  BLUE_SKY_TYPE_IMPL_SHORT(A, objbase, "A class")

	************************************************************************************

	PYTHON WRAPPER DEFINITION:

	Declarations:

	py_A : public py_named_pbase<A> {
		...
	};

	Export this with:

	py_exports_A () {
		py_export_named_pbase<A>("A_base");
    class_<py_A, bases< py_named_pbase<A> > >("A")
	}
*/

//! macro-definition for fast register element in property_base variant array
#define EREG(class_name,name,descr) ereg<class_name>(name, #name, descr)

#define EREG_OLOAD(class_name) \
	class_name &ereg(class_name ## _idxs idx, const char *name, const char *descr, named_pbase::etype it = named_pbase::PT_INT) { \
		this->params_names[idx] = named_pbase::ndt_t (std::string(name),std::string(descr),it); \
		this->assoc_names[name] = idx; \
		return *this; }

namespace blue_sky
  {
  /*!
   * \class named_pbase
   * \brief named-access wrapper for property_base
   */
  class BS_API_PLUGIN named_pbase : public property_base
    {
    public:
      //! indexing enumeration definition
      PROP_BASE_IDX_DECL_BEGIN(named_pbase,property_base)
      PROP_BASE_IDX_DECL_END

      //! access methods implementations
      PBASE_ACCESS_MS(named_pbase)
    public:
      enum etype
      {
        PT_NONE = -1,
        PT_INT,
        PT_FLOAT,
        PT_STR,
        PT_BOOL
      };

      struct ndt_t
        {
          ndt_t()
              : name(""), descr(""), type(PT_NONE)
          {}

          ndt_t(const ndt_t &src)
              : name(src.name), descr(src.descr), type(src.type)
          {}

          ndt_t(const std::string &n, const std::string &d, etype t)
              : name(n), descr(d), type(t)
          {}

          std::string name,   //!< name of property
          descr;  //!< description of property
          etype type;         //!< type of property
        };
      // typedefs
      typedef smart_ptr<named_pbase,true> sp_nprop; //!< smart_ptr of named_pbase
      //typedef std::pair<std::string, std::string> pair_t; //!< name and decription pair
      typedef std::vector< ndt_t > vndt_t; //!< type of name-descr-etype vector
      typedef std::map<std::string,int> container; //!< map for accessing by names instead indexes

      //! resize method
      void resize(int size)
      {
        property_base::resize(size);
        params_names.resize(size);
      }

      //! get name by index
      std::string get_name(idx_type idx) const
        {
          return params_names[idx].name;
        }

      //! get descrioption by index
      std::string get_descr(idx_type idx) const
        {
          return params_names[idx].descr;
        }

      etype get_etype(idx_type idx) const
        {
          return params_names[idx].type;
        }

      //! clear all arrays
      void clear()
      {
        assoc_names.clear();
        params_names.clear();
      }

      /*//! push_back method - may be necessary insert push_back::property_base();
      template<class T>
      void push_back(const pair_t &src) {
      	params_names::push_back(src);
      }*/

      const named_pbase &operator=(const named_pbase &src)
      {
        property_base::operator=(src);
        params_names = src.params_names;
        assoc_names = src.assoc_names;
        return *this;
      }

      const named_pbase &operator+=(const named_pbase &src)
      {
        property_base::operator+=(src);
        params_names.insert(params_names.end(),src.params_names.begin(),src.params_names.end());
        assoc_names.insert(src.assoc_names.begin(),src.assoc_names.end());
        return *this;
      }

      //! set value by name method
      template <class T>
      void 
      set_by_name (const std::string &name, const T &val)
      {
        container::iterator iter = assoc_names.find(name);
        if (iter == assoc_names.end())
          {
            bs_throw_exception ("No such element in map");
          }

        set_param <T> (iter->second, val, 0xdeadbeaf);
      }

      //! get value by name method
      template <class T>
      T 
      get_by_name(const std::string &name)  //, T t = T()) {
      {
        container::iterator iter = assoc_names.find(name);
        if (iter == assoc_names.end())
          {
            bs_throw_exception ("No such element in map");
          }

        return get_value_internal <T> (iter->second, 0xdeadbeaf);
      }

      bool
      set_value (int index, const std::string &value);

      size_t 
      size () const
      {
        return params_names.size ();
      }

      //! virtual dtor
      virtual ~named_pbase ();

      //! blue-sky necessary declarations
      BLUE_SKY_TYPE_DECL(named_pbase);

    protected:
      //! sets up names association with indexes
      void setup_assoc_names();

      //! register new element by index
      template<class T>
      named_pbase &ereg(typename T::idx_type idx, const char *name, const char *descr, etype it = PT_INT)
      {
        params_names[idx] = ndt_t (std::string(name),std::string(descr),it);
        assoc_names[name] = idx;
        return *this;
      }

      vndt_t params_names; //!< parameters names
    public:
      container assoc_names; //!< association map
    };

  namespace python
    {
#ifdef BSPY_EXPORTING_PLUGIN
    //! exports to pyhon. 
    inline void 
    py_export_named_pbase(const char *class_name) 
    {
			void (named_pbase::*set_by_name1)(const std::string&, const int &)				  = &named_pbase::set_by_name <int>;
      void (named_pbase::*set_by_name2)(const std::string&, const double&)	      = &named_pbase::set_by_name <double>;
      void (named_pbase::*set_by_name3)(const std::string&, const std::string&)   = &named_pbase::set_by_name <std::string>;
      void (named_pbase::*set_by_name4)(const std::string&, const bool&)				  = &named_pbase::set_by_name <bool>;

      using namespace boost::python;
      
      class_<named_pbase, boost::noncopyable> (class_name, no_init)
        .def ("__cons__",     make_constructor (construct_python_object <named_pbase>))
        .def ("__init__",     make_function (init_python_object <named_pbase>))
        .def("get_by_name_f", &named_pbase::get_by_name <double>)
        .def("get_by_name_i", &named_pbase::get_by_name <int>)
        .def("get_by_name_b", &named_pbase::get_by_name <bool>)
        .def("get_by_name_s", &named_pbase::get_by_name <std::string>)
        .def("set_by_name", set_by_name1)
        .def("set_by_name", set_by_name2)
        .def("set_by_name", set_by_name3)
        .def("set_by_name", set_by_name4)
        ;

      register_ptr_to_python <smart_ptr <named_pbase, true> > ();
    }


#endif
  } // end of namespace python
} // end of namespace blue_sky

#endif // NAMED_PBASE_ACCESS_H
