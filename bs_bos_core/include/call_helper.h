/**
 * \file   call_helper.cpp
 * \brief  Helper functions for call methods of derived from py_object classes
 * \author Miryanov Sergey
 * \date 2008-04-07
 */


#ifndef CALL_HELPER_H_
#define CALL_HELPER_H_

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky
  {
  namespace tx
    {

    using blue_sky::python::py_objbase;

    //////////////////////////////////////////////////////////////////////////
    template <class R, class C, class A>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A), A a)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a);
    }
    template <class R, class C, class A>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A&), A &a)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a);
    }

    //////////////////////////////////////////////////////////////////////////
    template <class R, class C, class A, class A2>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A, A2), A a, A2 a2)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a, a2);
    }
    template <class R, class C, class A, class A2>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A&, A2&), A &a, A2 &a2)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a, a2);
    }

    //////////////////////////////////////////////////////////////////////////
    template <class R, class C, class A, class A2, class A3>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A, A2, A3), A a, A2 a2, A3 a3)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a, a2, a3);
    }
    template <class R, class C, class A, class A2, class A3>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A&, A2&, A3&), A &a, A2 &a2, A3 &a3)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a, a2, a3);
    }
    template <class R, class C, class A, class A2, class A3, class A4>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A, A2, A3, A4), A a, A2 a2, A3 a3, A4 a4)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a, a2, a3, a4);
    }

    template <class R, class C, class A, class A2, class A3, class A4, class A5>
    R call(py_objbase * py, /*typename*/ R(C::*f)(A, A2, A3, A4), A a, A2 a2, A3 a3, A4 a4, A5 a5)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)(a, a2, a3, a4, a5);
    }

    template <class R, class C>
    R call (py_objbase * py, /*typename*/ R(C::*f)())
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)();
    }
    template <class R, class C>
    R call (py_objbase * py, /*typename*/ R(C::*f)() const)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)();
    }
    template <class R, class C>
    R call (const py_objbase * py, /*typename*/ R(C::*f)() const)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      return ((*c).*f)();
    }

    //////////////////////////////////////////////////////////////////////////
    template <class C, class A>
    void call(py_objbase * py, /*typename*/ void (C::*f)(A), A a)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      ((*c).*f)(a);
    }
    template <class C, class A>
    void call(py_objbase * py, /*typename*/ void (C::*f)(A&), A &a)
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      ((*c).*f)(a);
    }
    template <class C>
    void call(py_objbase * py, /*typename*/ void(C::*f)())
    {
      lsmart_ptr<smart_ptr<C, true> > c(py->get_sp());
      ((*c).*f)();
    }


  } // namespace tx
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef CALL_HELPER_H_
