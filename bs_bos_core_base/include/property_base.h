/**
 * @file property_base.h
 * @brief
 * @author Borschuk Oleg, modified by Nikonov Max
 * @date 2008-03-29
 */

#ifndef PROPERTY_BASE_H_
#define PROPERTY_BASE_H_

#include "bos_report.h"

//! calculate length
#define LENGTH_CALC(N,S) ((N) / (S) / 8 + 1)
//! find index of int variable in array for bit I
#define FIND_INDEX_IN_ULONG_ARRAY(I,S)                    ((int)((I) / (S) / 8))

//! find index in int variable for bit I
#define FIND_INDEX_IN_ULONG_VAR(I,S)                      ((I) - ((int)((I) / (S) / 8)) * (S) * 8)

//! return non zero if I-th bit in VAR is 1
#define CHECK_BIT_IN_ULONG_VAR(VAR,I)                   ((VAR) & (1 << (I)))

//! set to 1 the I-th bit in VAR
#define SET_BIT_IN_ULONG_VAR(VAR,I)                     ((VAR) = ((VAR) | (1 << (I))))

//! set to 0 the I-th bit in VAR
#define UNSET_BIT_IN_ULONG_VAR(VAR,I)                   ((VAR) = ((VAR) & (~(1 << (I)))))

//! return non zero if I-th bit in int array ARRAY is 1
#define CHECK_BIT_IN_ULONG_ARRAY(ARRAY,I,S) CHECK_BIT_IN_ULONG_VAR ((ARRAY)[FIND_INDEX_IN_ULONG_ARRAY (I,S)], FIND_INDEX_IN_ULONG_VAR (I,S))

//! set to 1 the I-th bit in int array ARRAY
#define SET_BIT_IN_ULONG_ARRAY(ARRAY,I,S)   SET_BIT_IN_ULONG_VAR ((ARRAY)[FIND_INDEX_IN_ULONG_ARRAY (I,S)], FIND_INDEX_IN_ULONG_VAR (I,S))

//! set to 0 the I-th bit in int array ARRAY
#define UNSET_BIT_IN_ULONG_ARRAY(ARRAY,I,S)   UNSET_BIT_IN_ULONG_VAR ((ARRAY)[FIND_INDEX_IN_ULONG_ARRAY (I,S)], FIND_INDEX_IN_ULONG_VAR (I,S))

namespace blue_sky
  {

  /*!
   * \brief boolvect ... short description ...
   * \author Oleg Borschuk <borschukos@ufantc.ru>
   * \date 2007-03-20
   * ... description ...
   * \comment boolvect doesn't work under "Linux geostation 2.6.27-10-generic #1 SMP Fri Nov 21 19:19:18 UTC 2008 x86_64 GNU/Linux)"
   *          this is a reason to wrap std::vector <bool> with bvector
   *          (Sergey Miryanov at 27.01.2009)
   */
    class bvector : protected std::vector <bool>
  {
  public:
    bvector () {}

    void 
    resize (int new_size)
    {
      assign (new_size, false);
    }

    //! clear buffer (unset all)
    void 
    clear ()
    {
      assign (size (), false);
    }

    //! return 0 if bit with index is unset
    bool 
    get (const int index) const
    {
      return this->operator[] (index);
    }

    //! return 0 if bit with index is unset
    void 
    set (const int index)
    {
      this->operator[] (index) = true;
    }
    //! return 0 if bit with index is unset
    void 
    unset (const int index)
    {
      this->operator[] (index) = false;
    }

    //bvector& operator= (const bvector &rhs)
    //{
    //  std::vector <bool>.operator = (rhs);
    //  return *this;
    //}
  };

  /*!
   * \namespace allowed_types
   * \brief Namespace that contains compile-time errors
   * rising procedures - for checking allowed in
   * property_base types.
   */
  namespace allowed_types
    {
    //! allowed types are int,float,bool,string
    BS_API_PLUGIN void cast_array_type(const int&);
    BS_API_PLUGIN void cast_array_type(const double&);
    BS_API_PLUGIN void cast_array_type(const bool&);
    BS_API_PLUGIN void cast_array_type(const std::string &);

    //! not allowed types are all other
    template<class T>
    void cast_array_type(const T&)
    {
      struct wrong_type;
      wrong_type smart_array_T;
    }

    /*template <class T, bool equ> struct equal_check;

    template<class T>
    struct equal_check<T, true> {
    	typedef T tval;
    };

    template<class T>
    struct equal_check<T, false> {
    	typedef T tval;
    	equal_check() {
    		struct incomplete_type;
    		incomplete_type smart_array_T;
    	}
    };*/
  }

// macro-definitions for making code writing simpler
//! begin of property indexing enumeration
//! @cname - inherited class name
//! @bcname - base class name - property_base for p-base-inherited types
#define PROP_BASE_IDX_DECL_BEGIN(cname,bcname) \
	public: enum cname ## _idxs { \
		__prop_base_idx_beg__ = (bcname::__prop_base_idx_end__ - 1),

//! end of property indexing enumeration
#define PROP_BASE_IDX_DECL_END \
		__prop_base_idx_end__ };

//! methods of property_base-inherited for access map of variants
// @class_name - name of property_base-inherited class
// get_param - returns value from container (idx - enum-based index)
// check_value - checks element of container (idx - enum-based index)
// set_param - sets element of container (idx enum-based index, value - allowed_type)
// operator[] - returns const reference to element in container (idx - enum-based index)
// get_ref - returns reference to element in container (idx - enum-based index)
#define PBASE_ACCESS_MS(class_name)																				\
	public:																																	\
    typedef double fp_type;                                               \
		typedef class_name ## _idxs idx_type;																	\
    inline int check_value (size_t idx, size_t /*dummy*/) const           \
    {													                                            \
			return this->is_inited.get ((int)idx);                              \
    }	                                                                    \
    template <class T> T get_value_internal (size_t idx, size_t /*dummy*/) const          \
    {                                                                     \
      BS_ASSERT (check_value (idx, 0xdeadbeaf)) (idx);                    \
      if ((int)idx >= this->get_length () || (int)idx < 0)                \
        {							                                                    \
				  BS_ASSERT (false && "index out of range");                      \
        }											                                            \
			allowed_types::cast_array_type(T());																\
      static T dummy__ = T ();                                            \
			return my_variant::get (this->values[idx], dummy__);                \
    }																			                                \
		template<class T> inline T get_param (idx_type idx) const             \
    {						                                                          \
			return get_value_internal <T> (idx, 0xdeadbeaf);                    \
    }												                                              \
		template<class T>	inline void set_param (idx_type idx, const T &value)\
    {			                                                                \
      my_variant::set (this->values[idx], value);											    \
			this->is_inited.set ((int)idx);                                     \
    }																		                                  \
    template<class T>	inline void set_param (size_t idx, const T &value, size_t /*dummy*/)\
    {			                                                                \
      my_variant::set (this->values[idx], value);											    \
			this->is_inited.set ((int)idx);                                     \
    }												                                              \
		inline int check_value (idx_type idx) const                           \
    {													                                            \
			return this->is_inited.get ((int)idx);                              \
    }														                                          \
    bool get_bool (idx_type idx) const {														      \
      BS_ASSERT (check_value (idx)) (idx);                                \
      return get_param <bool> (idx); }																		\
    int get_int (idx_type idx) const {															      \
      BS_ASSERT (check_value (idx)) (idx);                                \
      return get_param <int> (idx); }																		  \
    fp_type get_float (idx_type idx) const {													    \
      BS_ASSERT (check_value (idx)) (idx);                                \
      return get_param <fp_type> (idx); }																	\
    std::string get_str (idx_type idx) const {							              \
      BS_ASSERT (check_value (idx)) (idx);                                \
      return get_param <std::string> (idx); }														  \
    void set_bool (idx_type idx, bool v) {                                \
      set_param <bool> (idx, v); }                                        \
    void set_int (idx_type idx, int v) {                                  \
      set_param <int> (idx, v); }                                         \
    void set_float (idx_type idx, fp_type v) {                            \
      set_param <fp_type> (idx, v); }                                     \
    void set_str (idx_type idx, const std::string &v) {                   \
      set_param <std::string> (idx, v); }                                 \
    template <typename T>                                                 \
    void set_param_internal (idx_type idx, const T &t) {                  \
      set_param <T> (idx, t); }                                           \
    bool get_bool_d (idx_type idx, const bool &def = false) const         \
    { return check_value (idx) ? get_bool (idx) : def; }                  \
    int get_int_d (idx_type idx, const int &def = 0) const                \
    { return check_value (idx) ? get_int (idx) : def; }                   \
    fp_type get_float_d (idx_type idx, const fp_type &def = 0) const      \
    { return check_value (idx) ? get_float (idx) : def; }                 \
    std::string get_str_d (idx_type idx, const std::string &def = "") const \
    { return check_value (idx) ? get_str (idx) : def; }                   

      struct my_variant
      {
        typedef double fp_type;

        enum value_type
        {
          value_int, 
          value_float, 
          value_str, 
          value_bool, 
          value_none,
        };

        my_variant ()
        : type_ (value_none)
        , value_float_ (0)
        {

        }

        static int
        get (const my_variant &var, const int &unused)
        {
          BS_ASSERT (var.type_ == get_type (unused));
          return var.value_int_;
        }
        static fp_type
        get (const my_variant &var, const fp_type &unused)
        {
          BS_ASSERT (var.type_ == get_type (unused));
          return var.value_float_;
        }
        static const std::string 
        get (const my_variant &var, const std::string &unused)
        {
          BS_ASSERT (var.type_ == get_type (unused));
          return var.value_str_;
        }
        static bool 
        get (const my_variant &var, const bool &unused)
        {
          BS_ASSERT (var.type_ == get_type (unused));
          return var.value_bool_;
        }

        static void
        set (my_variant &var, const int &value)
        {
          if (var.type_ == value_none)
            var.type_ = get_type (value);

          BS_ASSERT (var.type_ == get_type (value));
          var.value_int_ = value;
        }
        static void
        set (my_variant &var, const fp_type &value)
        {
          if (var.type_ == value_none)
            var.type_ = get_type (value);

          BS_ASSERT (var.type_ == get_type (value));
          var.value_float_ = value;
        }
        static void
        set (my_variant &var, const std::string &value)
        {
          if (var.type_ == value_none)
            var.type_ = get_type (value);

          BS_ASSERT (var.type_ == get_type (value));
          var.value_str_ = value;
        }
        static void
        set (my_variant &var, const bool &value)
        {
          if (var.type_ == value_none)
            var.type_ = get_type (value);

          BS_ASSERT (var.type_ == get_type (value));
          var.value_bool_ = value;
        }

        static value_type
        get_type (const int &)
        {
          return value_int;
        }
        static value_type
        get_type (const fp_type &)
        {
          return value_float;
        }
        static value_type
        get_type (const std::string &)
        {
          return value_str;
        }
        static value_type
        get_type (const bool &)
        {
          return value_bool;
        }

        value_type    type_;
        union
        {
          int         value_int_;
          fp_type     value_float_;
          bool        value_bool_;
        };
        std::string   value_str_;
      };
  /*!
   * \class property_base
   * \brief base class of all properties
   */
  class BS_API_PLUGIN property_base: public objbase
    {
      //------------------------------------
      // METHODS
      //====================================

    protected:

    public:
      //! typedefs
      //typedef boost::variant<int,float,std::string,bool> var_t; // element of container
      typedef my_variant                      var_t;
      typedef std::vector< var_t >            vparams_t; // container
      typedef smart_ptr<property_base, true>  sp_pb; // sp

      //! enumeration of indexes
      enum property_base_idxs
      {
        __prop_base_idx_beg__ = 0,
        __prop_base_idx_end__
      };

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL(property_base);
      // Default destrructor
      virtual ~property_base ()
      {
        clear ();
      }

      //! access container elements methods
      PBASE_ACCESS_MS(property_base)

    public:
      //! get_params_name virtual method - will be overloaded in
      //! named_pbase_access (file: named_pbase_access.h)
      virtual const std::string &get_params_name (idx_type /*idx*/)
      {
        static std::string s("NULL");
        return s;
      }

      //! set default values
      virtual void set_default_values () {}

      //! resize methods
      void resize (int size)
      {
        values.resize (size);
        is_inited.resize (size);
        //clear ();
      }

      /*//! push back method - deprecated - will be removed in next versions
      template<class T>
      void push_back(const T &src) {
      	allowed_types::cast_array_type(T()); // test for allowed types
      	values.push_back(src);
      	//is_inited.push_back()
      }*/

      //! clear all
      void clear ()
      {
        values.clear();
        is_inited.clear ();
      }

      //! copy operator
      property_base& operator= (const property_base &rhs)
      {
        values = rhs.values;
        is_inited = rhs.is_inited;
        return *this;
      }
      //! update operator
      property_base& operator+= (const property_base &rhs)
      {
        BS_ASSERT ("In operator += " && (get_length() != rhs.get_length()) && "can be used for arrays of equil size only.");
        for (int i = 0; i < get_length(); ++i)
          {
            if (rhs.check_value ((property_base_idxs)i))
              {
                values[i] = rhs.values[i];
                is_inited.set (i);
              }
          }
        return *this;
      }

      //! unsets inited array
      void unset (const int idx)
      {
        is_inited.unset (idx);
      }

      //! return length
      int get_length () const
        {
          return (int)values.size ();
        }

    protected:
      bvector is_inited; //! inited check vector
      vparams_t values; //! vector of allowed types
    };

  typedef property_base::sp_pb sp_pb;
}//ns blue_sky



#endif // PROPERTY_BASE_H_
