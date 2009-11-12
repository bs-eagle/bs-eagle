/**
 *       \file  py_event_base.h
 *      \brief  Python wrappers for event_base, python events iterator,
 *              for event_base see event_base.h
 *     \author  Nikonov Max
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */
#ifndef PY_EVENT_BASE_H
#define PY_EVENT_BASE_H

#include "event_base.h"
#include "event_manager.h"

#include "calc_model.h"
#include "reservoir.h"

#include "export_python_wrapper.h"

#include "py_reservoir.h"
#include "py_calc_model.h"

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"

namespace blue_sky {
namespace python {

    STRATEGY_CLASS_WRAPPER (event_base, py_event_base)
    {
    public:
      typedef rs_mesh_iface <strategy_t>        mesh_t;
      typedef reservoir <strategy_t>            reservoir_t;
      typedef py_calc_model <strategy_t>        py_calc_model_t;
      typedef smart_ptr <mesh_t, true>          sp_mesh_iface_t;
      typedef smart_ptr <reservoir_t, true>     sp_reservoir_t;
      typedef event_base <strategy_t>           base_t;

    public:
      MAKE_ME_HAPPY (py_event_base, py_event_base_base <strategy_t>, "py_event_base");
      //WRAPPER_METHOD (apply, void, 3, (const sp_reservoir_t &, const py_mesh_iface_t &, const py_calc_model_t &));
    };

   // template <typename strategy_t>
   // class BS_API_PLUGIN py_event_base : public py_event_base_pv_holder<strategy_t>
      //                                , public boost::python::wrapper< py_event_base_pv_holder<strategy_t> > //objbase > 
   // {
   // public:
   //   typedef event_base <strategy_t>                 wrapped_t;
   //   //typedef py_named_pbase                          base_t;
   //   typedef smart_ptr < wrapped_t >                 sp_event_base_t;
   //   typedef py_event_base <strategy_t>              this_t;

   //   typedef py_reservoir <strategy_t>               py_reservoir_t;
   //   typedef py_rs_mesh <strategy_t>                 py_mesh_iface_t;
   //   typedef py_calc_model <strategy_t>              py_calc_model_t;

   //   typedef smart_ptr < reservoir <strategy_t> >    sp_top_t;
   //   typedef smart_ptr < rs_mesh_iface <strategy_t> >      sp_mesh_iface_t;
   //   typedef smart_ptr < calc_model <strategy_t> >   sp_calc_model_t;

   //   py_event_base ();
   //   py_event_base (const sp_obj &);
   //   py_event_base (const py_objbase &);

   //   void init(const std::string & _params);

   //   //virtual void parse();
   //   void py_apply (const py_reservoir_t &top,
   //                  const py_mesh_iface_t &mesh,
   //                  const py_calc_model_t &calc_model);

   //   void apply (const sp_top_t &top,
   //               const sp_mesh_iface_t &mesh,
   //               const sp_calc_model_t &calc_model);

   //   void apply_from_py (const py_reservoir_t &top,
   //                       const py_mesh_iface_t &mesh,
   //                       const py_calc_model_t &calc_model);

      //void default_apply_from_py (const py_reservoir_t &top,
      //                            const py_mesh_iface_t &mesh,
      //                            const py_calc_model_t &calc_model);
   // };


    template <typename strategy_t> class py_el_pair;

    template <typename strategy_t>
    class py_event_base_iterator : public std::iterator<
          std::bidirectional_iterator_tag,
          py_event_base<strategy_t>, ptrdiff_t,
          py_event_base<strategy_t>, py_event_base<strategy_t> >
    {

    public:

      typedef event_manager <strategy_t>           event_manager_t;
      typedef typename event_manager_t::event_map  em_t;
      typedef typename em_t::iterator              list_iterator_t;
      typedef typename list_iterator_t::value_type value_t;
      typedef typename value_t::second_type        elist_t;

      typedef py_el_pair <strategy_t>              py_el_pair_t;
      typedef typename elist_t::value_type         sp_event_base_t;
      typedef typename elist_t::iterator           iterator_t;
      typedef typename elist_t::const_iterator     const_iterator_t;

      typedef py_event_base <strategy_t>           py_event_base_t;
      typedef py_event_base_iterator <strategy_t>  this_t;

      typedef std::iterator<std::bidirectional_iterator_tag,
      py_event_base<strategy_t>, ptrdiff_t,
      py_event_base<strategy_t>, py_event_base<strategy_t> >  base_t;

      typedef typename base_t::reference  reference;
      typedef typename base_t::pointer    pointer;

      py_event_base_iterator (const elist_t &events);
      py_event_base_iterator (const this_t &iter);

      py_event_base_t *next ();

      reference operator*() const;

      this_t& operator++();
      this_t operator++(int);

      this_t& operator--();
      this_t operator--(int);

      bool operator ==(const this_t &ritr) const;
      bool operator !=(const this_t &ritr) const;
      const this_t &operator =(const this_t &ritr);

    protected:
      elist_t events;
      const_iterator_t iter;
    };

    template <typename strategy_t>
    class py_el_pair
    {
    public:
      typedef py_el_pair <strategy_t>       this_t;

      typedef event_manager <strategy_t>           event_manager_t;
      typedef typename event_manager_t::event_map  em_t;
      typedef typename em_t::const_iterator        const_iterator_t;
      typedef typename em_t::iterator              iterator_t;
      typedef typename iterator_t::value_type      value_t;
      typedef typename value_t::first_type         first_t; // posix_time::ptime
      typedef typename value_t::second_type        second_t; // list < smart_ptr <event_base> >

      typedef py_event_base_iterator <strategy_t>  py_event_base_iterator_t;

      py_el_pair ();
      py_el_pair (const const_iterator_t &iter);
      py_el_pair (const this_t &iter);

      py_event_base_iterator_t list_begin ();

    protected:
      first_t  first;
      second_t second;
      //iterator_t &iter;
    };

    template <typename strategy_t>
    class event_list_iterator : public std::iterator<
          std::bidirectional_iterator_tag,
          py_el_pair<strategy_t>, ptrdiff_t,
          py_el_pair<strategy_t>, py_el_pair<strategy_t> >
    {

    public:
      typedef event_list_iterator <strategy_t>    this_t;

      typedef event_manager <strategy_t>           event_manager_t;
      typedef smart_ptr <event_manager_t>          sp_event_manager_t;
      typedef typename event_manager_t::event_map  em_t;
      typedef typename em_t::const_iterator        const_iterator_t;
      typedef typename em_t::iterator              iterator_t;
      typedef typename em_t::key_type              first_t;
      typedef typename em_t::value_type            second_t;
      typedef py_el_pair <strategy_t>              py_el_pair_t;

      typedef std::iterator<std::bidirectional_iterator_tag,
      py_el_pair<strategy_t>, ptrdiff_t,
      py_el_pair<strategy_t>, py_el_pair<strategy_t> >  base_t;

      typedef typename base_t::reference  reference;
      typedef typename base_t::pointer    pointer;


      event_list_iterator (const sp_event_manager_t &); //iterator_t &iter);
      event_list_iterator (const this_t &iter);

      reference operator*() const;

      this_t& operator++();
      this_t operator++(int);

      py_el_pair_t next ();

      this_t& operator--();
      this_t operator--(int);

      bool operator ==(const this_t &ritr) const;
      bool operator !=(const this_t &ritr) const;
      const this_t &operator =(const this_t &ritr);

    protected:
      sp_event_manager_t evm;
      const_iterator_t  iter;
      //py_el_pair_t pair_;
    };

    /**
     * \brief  Exports wrappers to python
     * */
    void 
    py_export_events ();

} // namespace python
} // namespace blue_sky

#endif // PY_EVENT_BASE_H
