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
   */
  class bvector
    {
      typedef unsigned long ulong;
      typedef std::vector<ulong> ulong_v;
      //-----------------------------------
      // METHODS
      //-----------------------------------
    public:
      // Default constructor
      bvector () {}

      // Destructor
      virtual ~bvector ()
      {
        clear ();
      }

      void resize (const int size)
      {
        buffer_.resize (LENGTH_CALC (size, sizeof (ulong)));
        clear ();
      }

      //! return const reference to the buffer
      const ulong_v &buffer () const
        {
          return buffer_;
        }

      //! return buffer length in unsigned long
      int buffer_length_in_ulong ()
      {
        return (int)buffer_.size ();
      }

      //! return buffer length in bytes
      int buffer_length_in_bytes ()
      {
        return (int)buffer_.size () * sizeof (ulong);
      }

      //! clear buffer (unset all)
      void clear ()
      {
        memset (&buffer_[0], 0, sizeof (ulong) * buffer_.size ());
      }

      //! return 0 if bit with index is unset
      int get (const int index) const
        {
          return (CHECK_BIT_IN_ULONG_ARRAY (buffer_, index, sizeof (ulong)));
        }

      //! return 0 if bit with index is unset
      void set (const int index)
      {
        (SET_BIT_IN_ULONG_ARRAY (buffer_, index, sizeof (ulong)));
      }
      //! return 0 if bit with index is unset
      void unset (const int index)
      {
        (UNSET_BIT_IN_ULONG_ARRAY (buffer_, index, sizeof (ulong)));
      }

      bvector& operator= (const bvector &rhs)
      {
        buffer_ = rhs.buffer_;
        return *this;
      }

    protected:
    private:
      //-----------------------------------
      // VARIABLES
      //-----------------------------------
    public:
    protected:
    private:
      std::vector<ulong> buffer_;
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
    BS_API_PLUGIN void cast_array_type(const float&);
    BS_API_PLUGIN void cast_array_type(const bool&);
    BS_API_PLUGIN void cast_array_type(const std::string&);

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
#define PB_IDX_BEG(cname,bcname,type) \
	public: enum cname ## _idxs_ ## type { \
		pb_idx_ ## type ## _beg = (bcname::pb_idx_ ## type ## _end - 1),

//! end of property indexing enumeration
#define PB_IDX_END(type) \
		pb_idx_ ## type ## _end };

//! methods of property_base-inherited for access map of variants
// @class_name - name of property_base-inherited class
// get_param - returns value from container (idx - enum-based index)
// check_value - checks element of container (idx - enum-based index)
// set_param - sets element of container (idx enum-based index, value - allowed_type)
// operator[] - returns const reference to element in container (idx - enum-based index)
// get_ref - returns reference to element in container (idx - enum-based index)
#define PB_IMPL(class_name) \
	public: \
		typedef class_name ## _idxs_int idx_i_type; \
		typedef class_name ## _idxs_float idx_f_type;	\
		typedef class_name ## _idxs_string idx_s_type;	\
		typedef class_name ## _idxs_bool idx_b_type; \
		const int &get_param(idx_i_type idx) const { return this-> /*template*/ get_param(idx,this->ivals); } \
		const float &get_param(idx_f_type idx) const { return this-> /*template*/ get_param(idx,this->fvals); } \
		const char &get_param(idx_b_type idx) const { return this-> /*template*/ get_param(idx,this->bvals); } \
		const std::string &get_param(idx_s_type idx) const { return this-> /*template*/ get_param(idx,this->svals); } \
		\
		int &get_param(idx_i_type idx) { return this-> /*template*/ get_param(idx,this->ivals); } \
		float &get_param(idx_f_type idx) { return this-> /*template*/ get_param(idx,this->fvals); } \
		char &get_param(idx_b_type idx) { return this-> /*template*/ get_param(idx,this->bvals); } \
		std::string &get_param(idx_s_type idx) { return this-> /*template*/ get_param(idx,this->svals); } \
		\
		inline int check_value (idx_i_type idx) const { return this-> /*template*/ check_value(idx,this->is_inited_i); } \
		inline int check_value (idx_f_type idx) const { return this-> /*template*/ check_value(idx,this->is_inited_f); } \
		inline int check_value (idx_b_type idx) const { return this-> /*template*/ check_value(idx,this->is_inited_b); } \
		inline int check_value (idx_s_type idx) const { return this-> /*template*/ check_value(idx,this->is_inited_s); } \
		\
		inline void set_param (idx_i_type idx, const int &value) { this-> /*template*/ set_param(idx,this->ivals,value,this->is_inited_i); } \
		inline void set_param (idx_f_type idx, const float &value) { this-> /*template*/ set_param(idx,this->fvals,value,this->is_inited_f); } \
		inline void set_param (idx_b_type idx, const bool &value) { this-> /*template*/ set_param(idx,this->bvals,(char)value,this->is_inited_b); } \
		inline void set_param (idx_s_type idx, const std::string &value) { this-> /*template*/ set_param(idx,this->svals,value,this->is_inited_s); } \
		\
		inline void unset (idx_i_type idx) { this-> /*template*/ unset(idx,this->is_inited_i); } \
		inline void unset (idx_f_type idx) { this-> /*template*/ unset(idx,this->is_inited_f); } \
		inline void unset (idx_b_type idx) { this-> /*template*/ unset(idx,this->is_inited_b); } \
		inline void unset (idx_s_type idx) { this-> /*template*/ unset(idx,this->is_inited_s); }

  /*!
   * \class property_base
   * \brief base class of all properties
   */
  class BS_API_PLUGIN property_base: public objbase
    {
      //------------------------------------
      // METHODS
      //====================================
    public:
      //! typedefs
      //typedef boost::variant<int,float,std::string,bool> var_t; // element of container
      //typedef std::vector< var_t > vparams_t; // container
      typedef std::vector< int >             vparams_i; //!< int params
      typedef std::vector< float >           vparams_f; //!< float params
      typedef std::vector< char >            vparams_b; //!< bool params
      typedef std::vector< std::string >     vparams_s; //!< string params
      typedef smart_ptr<property_base, true> sp_pb;     //!< sp

      //! enumeration of int indexes
      enum property_base_idxs_int
      {
        pb_idx_int_beg = 0,
        pb_idx_int_end
      };

      //! enumeration of float indexes
      enum property_base_idxs_float
      {
        pb_idx_float_beg = 0,
        pb_idx_float_end
      };

      //! enumeration of string indexes
      enum property_base_idxs_string
      {
        pb_idx_string_beg = 0,
        pb_idx_string_end
      };

      //! enumeration of bool indexes
      enum property_base_idxs_bool
      {
        pb_idx_bool_beg = 0,
        pb_idx_bool_end
      };

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL(property_base);
      // Default destrructor
      virtual ~property_base ()
      {
        clear ();
      }

      //! access container elements methods
      PB_IMPL(property_base)

    public:
      //! get_params_name virtual method - will be overloaded in
      //! named_pbase_access (file: named_pbase_access.h)
      /*virtual const std::string &get_params_name (idx_type) {
      	static std::string s("NULL");
      	return s;
      }*/

      //! set default values
      virtual void set_default_values () {}

      //! resize methods
      void resize_int (int size);
      void resize_float (int size);
      void resize_bool (int size);
      void resize_string (int size);

      //! clear all
      void clear ();

      //! copy operator
      property_base& operator= (const property_base &rhs);
      //! update operator
      property_base& operator+= (const property_base &rhs);

    protected:
      template <class T, class R>
      const typename R::value_type &get_param(T idx, const R &r) const
        {
          BS_ASSERT((int)r.size() > (int)idx);
          return r[idx];
        }

      template <class T, class R>
      typename R::value_type &get_param(T idx, R &r)
      {
        BS_ASSERT((int)r.size() > (int)idx);
        return r[idx];
      }

      template<class T, class R>
      const R& operator[](T idx) const
        {
          return this-> /*template*/ get_param(idx);
        }

      template<class T, class R>
      R& operator[](T idx)
      {
        return this-> /*template*/ get_param(idx);
      }

      template <class T>
      inline int check_value(T idx, const bvector &v) const
        {
          return v.get((int)idx);
        }

      template <class T>
      inline void unset(T idx, bvector &v)
      {
        v.unset((int)idx);
      }

      template <class V, class I>
      void resize (int size, V &v, I &i)
      {
        v.resize (size);
        i.resize (size);
      }

      template<class T, class R>
      inline void set_param (T idx, R &params, const typename R::value_type &value, bvector &is_inited)
      {
        params[idx] = value;
        is_inited.set ((int)idx);
      }

      bvector is_inited_i, //! inited check vector int
      is_inited_f, //! inited check vector float
      is_inited_b, //! inited check vector bool
      is_inited_s; //! inited check vector string

      //vparams_t values; //! vector of allowed types
      vparams_i ivals; //!< vector of integers
      vparams_f fvals; //!< vector of floats
      vparams_b bvals; //!< vector of bools
      vparams_s svals; //!< vector of strings
    };

  typedef property_base::sp_pb sp_pb;
}//ns blue_sky



#endif // PROPERTY_BASE_H_
