#ifndef ARITHMETIC_OPERAT_H
#define ARITHMETIC_OPERAT_H

#include "ar_stack.h"

namespace blue_sky
  {

	template<typename strategy_t>
  class BS_API_PLUGIN ar_operat
    {
    public:
      typedef ar_args <strategy_t>   ar_args_t;
			typedef ar_stack <strategy_t>  ar_stack_t;
			typedef ar_operat <strategy_t> this_t;

			typedef typename strategy_t::item_t item_t;

      typedef int (*operat_func) (ar_args_t &ar,  ar_stack <ar_args_t> &st);

      ar_operat();
      ar_operat(const std::string &tname);
      ar_operat(operat_func toper, int tpriority, int tnum_of_arg, const std::string &tname);

      virtual ~ar_operat();

      void init();

      int set_name(const std::string &tname);
      const std::string &get_name() const;

      int operator==(const std::string &tname) const;
      bool operator<(const ar_operat &op) const;
      //int operator() (ar_args<T> &ar, std::stack<T> &st) = 0;

      //! data
      operat_func oper;
      int priority,num_of_arg;

    private:
      std::string name;
    };

#define A_GET_VALUE(a,n)     \
	((a.get_type () == ar_args<strategy_t>::FPOINT_T) ? (((typename strategy_t::item_t *)(a.get_array ()))[n]) \
																				 : (((int   *)(a.get_array ()))[n]))

#define A_SET_VALUE(l,n,value)  \
	((l.get_type () == ar_args<strategy_t>::FPOINT_T) ?  (((typename strategy_t::item_t *)(l.get_array ()))[n] = (typename strategy_t::item_t)(value)) \
																				 :  (((int   *)(l.get_array ()))[n] = (int)(value)))

  // operators and function

  /*!
  	\brief check function for many operations
  */

  /*!
  	\brief all_operations
  */
  //template<typename T>
	template <typename strategy_t, typename T>
	int all_trigon_operations(T(*m_func)(T),
                            ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st, bool check = false,
														typename strategy_t::item_t val = 0)
  {
    ar_args<strategy_t> l = st.pop();
    if (l.def_flag)
      {
        std::ostringstream msg;
        msg << "Value of token " << l.get_name() << " is undefined";
        throw bs_exception("ar_operation",msg.str().c_str());
      }
    ar_args<strategy_t> res(l);

    if (res.ni * res.nj * res.nk == 1)
      {
        A_SET_VALUE (res,0,m_func (A_GET_VALUE(res,0)));
      }
    else
      {
        int index = -1;
        for (int k = 0; k < res.nk; ++k)
          for (int j = 0; j < res.nj; ++j)
            for (int i = 0; i < res.ni; ++i)
              {
                index = i + j*res.ni + k*res.ni*res.nj;
                if (check)
                  //if (res[index] == val) {
                  if (A_GET_VALUE(res,index) == val)
                    {
                      std::ostringstream msg;
                      msg << "argument shouldn't be equal " << val << "!";
                      throw bs_exception("operation",msg.str().c_str());
                      //return YS_BAD_VALUE;
                    }
                //res[index] = m_func(res[index]);
                A_SET_VALUE (res,index,m_func (A_GET_VALUE(res,index)));
              }
      }
    st.push(res);
    return 0;
  }

  //template<typename T>
	template <typename strategy_t>
  int ar_operation_abs(ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_trigon_operations <strategy_t, double> (fabs,ar,st);
  }

  //template<typename T>
	template <typename strategy_t>
  int ar_operation_exp (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_trigon_operations <strategy_t, double> (exp,ar,st);
  }

	//template<typename T>
	template <typename strategy_t>
  int ar_operation_cos (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_trigon_operations <strategy_t, double> (cos,ar,st);
  }

	//template<typename T>
	template <typename strategy_t>
  int ar_operation_sin (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_trigon_operations <strategy_t, double> (sin,ar,st);
  }

	//template<typename T>
	template <typename strategy_t>
  int ar_operation_log (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_trigon_operations <strategy_t, double> (::log,ar,st);
  }

	//template<typename T>
	template <typename strategy_t>
  int ar_operation_log10 (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_trigon_operations <strategy_t, double> (log10,ar,st);
  }

	//template<typename T>
	template <typename strategy_t>
  int ar_operation_sqrt (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_trigon_operations <strategy_t, double> (sqrt,ar,st);
  }

	//template<typename T>
	template <typename strategy_t>
  int ar_operation_tan (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
		return all_trigon_operations <strategy_t, double> (tan,ar,st,true,(((typename strategy_t::item_t)3.141592654) / 2));
  }

  //! math operations
  template<typename T>
  T add_operation (const T &l, const T &r)
  {
    return (l + r);
  }

  template<typename T>
  T sub_operation (const T &l, const T &r)
  {
    return (l - r);
  }

  template<typename T>
  T mul_operation (const T &l, const T &r)
  {
    return (l * r);
  }

  template<typename T>
  T dev_operation (const T &l, const T &r)
  {
    return (l / r);
  }

  template <typename strategy_t>
	int all_simple_operations (typename strategy_t::item_t(*math_operation)(const typename strategy_t::item_t&,const typename strategy_t::item_t&),
                             ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st,
														 bool check = false, typename strategy_t::item_t val = 0)
  {
    int fl_l = 0, fl_r = 0, index;
    ar_args<strategy_t> r = st.pop();
    ar_args<strategy_t> l = st.pop();
    if (l.def_flag)
      {
        std::ostringstream msg;
        msg << "Value of token " << l.get_name() << " is undefined";
        throw bs_exception("ar_operation",msg.str().c_str());
      }
    if (r.def_flag)
      {
        std::ostringstream msg;
        msg << "Value of token " << r.get_name() << " is undefined";
        throw bs_exception("ar_operation",msg.str().c_str());
      }
    ar_args<strategy_t> res(ar);

    if (r.ni * r.nj * r.nk == 1)
      fl_r = 1;
    else
      if (res.i2 > r.ni || res.j2 > r.nj || res.k2 > r.nk)
        throw bs_exception("ar_operation","Error: index of right side argument is out of range");
		if (l.ni * l.nj * l.nk == 1)
			fl_l = 1;
		else
			if (res.i2 > l.ni || res.j2 > l.nj || res.k2 > l.nk)
				throw bs_exception("ar_operation","Error: index of left side argument is out of range");

    if (fl_r == 1 && fl_l == 1)
      {
        res.ni = res.nj = res.nk = 1;
				if (check) {
					if (A_GET_VALUE(r,0) != val)
            A_SET_VALUE(res,0,math_operation(A_GET_VALUE(l,0),A_GET_VALUE(r,0)));
					else
						BOSERR (section::arithmetic, level::error) << "ar_operation internal error!" << bs_end;
				}	else {
					A_SET_VALUE(res,0,math_operation(A_GET_VALUE(l,0),A_GET_VALUE(r,0)));
				}
      }
    else if (fl_r == 1)
      {
        for (int k = res.k1; k < res.k2; ++k)
          for (int j = res.j1; j < res.j2; ++j)
            for (int i = res.i1; i < res.i2; ++i)
              {
								if (check) {
									if (A_GET_VALUE(r,0) != val)
                    {
                      index = i + j*res.ni + k*res.ni*res.nj;
                      A_SET_VALUE(res,index,math_operation(A_GET_VALUE(l,index),A_GET_VALUE(r,0)));
                    }
								}
								else {
									index = i + j*res.ni + k*res.ni*res.nj;
									A_SET_VALUE(res,index,math_operation(A_GET_VALUE(l,index),A_GET_VALUE(r,0)));
								}

              }
      }
    else if (fl_l == 1)
      {
        for (int k = res.k1; k < res.k2; ++k)
          for (int j = res.j1; j < res.j2; ++j)
            for (int i = res.i1; i < res.i2; ++i)
              {
                index = i + j*res.ni + k*res.ni*res.nj;
								if (check) {
									if (A_GET_VALUE(r,index) != val)
                    A_SET_VALUE(res,index,math_operation(A_GET_VALUE(l,0),A_GET_VALUE(r,index)));
								}
								else
									A_SET_VALUE(res,index,math_operation(A_GET_VALUE(l,0),A_GET_VALUE(r,index)));
              }
      }
    else
      {
        for (int k = res.k1; k < res.k2; ++k)
          for (int j = res.j1; j < res.j2; ++j)
            for (int i = res.i1; i < res.i2; ++i)
              {
                index = i + j*res.ni + k*res.ni*res.nj;
								if (check) {
									if (A_GET_VALUE(r,index) != val)
                    A_SET_VALUE(res,index,math_operation(A_GET_VALUE(l,index),A_GET_VALUE(r,index)));
								}
								else
									A_SET_VALUE(res,index,math_operation(A_GET_VALUE(l,index),A_GET_VALUE(r,index)));
              }
      }

		st.push (res);

    return 0;
  }

  template<typename strategy_t>
  int ar_operation_plus (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
		return all_simple_operations <strategy_t> (add_operation<typename strategy_t::item_t>,ar,st);
  }

  template<typename strategy_t>
  int ar_operation_minus (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_simple_operations <strategy_t> (sub_operation<typename strategy_t::item_t>,ar,st);
  }

  template<typename strategy_t>
  int ar_operation_mult (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_simple_operations <strategy_t> (mul_operation<typename strategy_t::item_t>,ar,st);
  }

  template<typename strategy_t>
  int ar_operation_dev (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_simple_operations <strategy_t> (dev_operation<typename strategy_t::item_t>,ar,st,true,0);
  }
  //! operation unary minus
  template<typename strategy_t>
  int ar_operation_rminus (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    int index;
    ar_args<strategy_t> l = st.pop();
    if (l.def_flag)
      {
        std::ostringstream msg;
        msg << "Value of token " << l.get_name() << " is undefined";
        throw bs_exception("ar_operation",msg.str().c_str());
      }
    ar_args<strategy_t> res(l);
    if (res.ni * res.nj * res.nk == 1)
      A_SET_VALUE(res,0,-A_GET_VALUE(res,0));
    else
      {
        for (int k = 0; k < res.nk; ++k)
          for (int j = 0; j < res.nj; ++j)
            for (int i = 0; i < res.ni; ++i)
              {
                index = i + j*res.ni + k*res.ni*res.nj;
                A_SET_VALUE(res,index,-A_GET_VALUE(res,index));
              }
      }
    return 0;
  }

#define GET_MAX(A,B) ((A)>(B)?(A):(B))
#define GET_MIN(A,B) ((A)<(B)?(A):(B))

  template<typename T>
  T operation_max (const T &l, const T &r)
  {
    return GET_MAX(l,r);
  }

  template<typename T>
  T operation_min (const T &l, const T &r)
  {
    return GET_MIN(l,r);
  }

  template<typename strategy_t>
  int ar_operation_max (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
		return all_simple_operations <strategy_t> (operation_max<typename strategy_t::item_t>,ar,st);
  }

  template<typename strategy_t>
  int ar_operation_min (ar_args<strategy_t> &ar, ar_stack < ar_args<strategy_t> > &st)
  {
    return all_simple_operations <strategy_t> (operation_min<typename strategy_t::item_t>,ar,st);
  }

  /*template <typename T>
  int all_compare_operations (T (*cmp_operation)(const T&,const T&), ar_args<T> &ar, ar_stack<T> &st) {
  	int fl_r = 0, fl_l = 0, index;
  	ar_args<T> r = st.pop();
  	ar_args<T> l = st.pop();
  	if (l.def_flag) {
  		std::ostringstream msg;
  		msg << "Value of token " << l.get_name() << " is undefined";
  		throw bs_exception("ar_operation",msg.str().c_str());
  	}
  	if (r.def_flag) {
  		std::ostringstream msg;
  		msg << "Value of token " << r.get_name() << " is undefined";
  		throw bs_exception("ar_operation",msg.str().c_str());
  	}
  	ar_args<T> res(ar);

  }*/
}

#endif // ARITHMETIC_OPERAT_H
