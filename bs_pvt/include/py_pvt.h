/**
 * \file py_pvt.h
 * \brief python wrappers for pvt_*
 * \author Miryanov Sergey
 * \date 08.05.2008
 */
#ifndef PY_BS_PVT_H_
#define PY_BS_PVT_H_

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky
  {
  namespace python
    {
    /**
     * \brief base class for python wrappers of pvt classes
     */
    class py_pvt_base : public py_objbase
    {
    public:

      /**
       * \brief constructor
       */
      py_pvt_base (sp_obj sp_obj_)
      : py_objbase (sp_obj_)
      {

      }

    };

    template <class pvt_t>
    class py_pvt_wrapper : public py_pvt_base
    {
    public:

      typedef py_pvt_base                            base_t;
      typedef py_pvt_wrapper<pvt_t>                  this_t;

      typedef typename pvt_t::pvt_strategy_t         pvt_strategy_t;
      typedef typename pvt_strategy_t::item_t        item_t;
      typedef typename pvt_strategy_t::item_array_t  item_array_t;

      typedef pvt_t                                  wrapped_t;

    public:

      /**
       * \brief constructor
       */
      py_pvt_wrapper ()
      : py_pvt_base (BS_KERNEL.create_object (pvt_t::bs_type ()))
      {

      }
      py_pvt_wrapper (const sp_obj &sp_obj_)
      : py_pvt_base (sp_obj_)
      {

      }

      void build (item_t atm_p, item_t min_p, item_t max_p, int n_intervals)
      {
        get_spx (this)->build (atm_p, min_p, max_p, n_intervals);
      }

      item_t get_p_step ()
      {
        return get_spx (this)->get_p_step ();
      }

      item_t get_surface_density ()
      {
        return get_spx (this)->get_surface_density ();
      }

      void set_surface_density (item_t surface_density)
      {
        get_spx (this)->set_surface_density (surface_density);
      }

    };

    template <typename pvt_t>
    class py_pvt_water : public py_pvt_wrapper <pvt_t> 
    {
    public:
      typedef py_pvt_wrapper <pvt_t>          base_t;
      typedef py_pvt_water <pvt_t>            this_t;
      typedef typename pvt_t::index_t         index_t;
      typedef typename pvt_t::item_t          item_t;
      typedef typename pvt_t::index_array_t   index_array_t;
      typedef typename pvt_t::item_array_t    item_array_t;
      typedef pvt_t                           wrapped_t;

    public:

      py_pvt_water (const sp_obj &sp)
      : base_t (sp)
      {
      }

      void
      calc (index_t n, index_t n_phases, item_t p, item_array_t &inv_fvf, item_array_t &p_inv_fvf,
        item_array_t &inv_visc, item_array_t &p_inv_visc,
        item_array_t &inv_visc_fvf, item_array_t &p_inv_visc_fvf)
      {
        BS_ASSERT (inv_fvf.size () && p_inv_fvf.size () && inv_visc.size () && p_inv_visc.size () && inv_visc_fvf.size () && p_inv_visc_fvf.size ()) (inv_fvf.size ()) (p_inv_fvf.size ()) (inv_visc.size ()) (p_inv_visc.size ()) (inv_visc_fvf.size ()) (p_inv_visc_fvf.size ());
        if (!inv_fvf.size () || !p_inv_fvf.size () || !inv_visc.size () || !p_inv_visc.size () || !inv_visc_fvf.size () || !p_inv_visc_fvf.size ())
          {
            throw bs_exception ("py_pvt_water::calc", "invalid argument");
          }

        index_t i = n * n_phases + 0;
        get_spx (this)->calc (p, &inv_fvf [i], &p_inv_fvf[i], &inv_visc[i], &p_inv_visc[i], &inv_visc_fvf[i], &p_inv_visc_fvf[i]);
      }
    };

    template <typename pvt_t>
    class py_pvt_gas : public py_pvt_wrapper <pvt_t> 
    {
    public:
      typedef py_pvt_wrapper <pvt_t>          base_t;
      typedef py_pvt_gas <pvt_t>              this_t;
      typedef typename pvt_t::index_t         index_t;
      typedef typename pvt_t::item_t          item_t;
      typedef typename pvt_t::index_array_t   index_array_t;
      typedef typename pvt_t::item_array_t    item_array_t;
      typedef pvt_t                           wrapped_t;

    public:

      py_pvt_gas (const sp_obj &sp)
      : base_t (sp)
      {
      }

      void
      calc (index_t n, index_t n_phases, item_t p, item_array_t &inv_fvf, item_array_t &p_inv_fvf,
        item_array_t &inv_visc, item_array_t &p_inv_visc,
        item_array_t &inv_visc_fvf, item_array_t &p_inv_visc_fvf)
      {
        BS_ASSERT (inv_fvf.size () && p_inv_fvf.size () && inv_visc.size () && p_inv_visc.size () && inv_visc_fvf.size () && p_inv_visc_fvf.size ()) (inv_fvf.size ()) (p_inv_fvf.size ()) (inv_visc.size ()) (p_inv_visc.size ()) (inv_visc_fvf.size ()) (p_inv_visc_fvf.size ());
        if (!inv_fvf.size () || !p_inv_fvf.size () || !inv_visc.size () || !p_inv_visc.size () || !inv_visc_fvf.size () || !p_inv_visc_fvf.size ())
          {
            throw bs_exception ("py_pvt_gas::calc", "invalid argument");
          }

        // TODO: BUG
        index_t i = n * n_phases + 1;
        get_spx (this)->calc (p, &inv_fvf [i], &p_inv_fvf[i], &inv_visc[i], &p_inv_visc[i], &inv_visc_fvf[i], &p_inv_visc_fvf[i]);
      }
    };

    template <typename pvt_t>
    class py_pvt_oil : public py_pvt_wrapper <pvt_t> 
    {
    public:
      typedef py_pvt_wrapper <pvt_t>          base_t;
      typedef py_pvt_oil <pvt_t>              this_t;
      typedef typename pvt_t::index_t         index_t;
      typedef typename pvt_t::item_t          item_t;
      typedef typename pvt_t::index_array_t   index_array_t;
      typedef typename pvt_t::item_array_t    item_array_t;
      typedef pvt_t                           wrapped_t;

    public:

      py_pvt_oil (const sp_obj &sp)
      : base_t (sp)
      {
      }

      void
      calc (index_t n, index_t n_phases, bool is_g, int main_var, item_t p, item_t gor,
        item_array_t &inv_fvf, item_array_t &p_inv_fvf, item_array_t &gor_inv_fvf,
        item_array_t &inv_visc, item_array_t &p_inv_visc, item_array_t &gor_inv_visc,
        item_array_t &inv_visc_fvf, item_array_t &p_inv_visc_fvf, item_array_t &gor_inv_visc_fvf,
        item_array_t &gas_oil_ratio, item_array_t &d_gas_oil_ratio, item_t drsdt, item_t dt,
        item_array_t &old_gas_oil_ratio)
      {
        BS_ASSERT (inv_fvf.size () && p_inv_fvf.size () && inv_visc.size () && p_inv_visc.size () && inv_visc_fvf.size () && p_inv_visc_fvf.size ()) (inv_fvf.size ()) (p_inv_fvf.size ()) (inv_visc.size ()) (p_inv_visc.size ()) (inv_visc_fvf.size ()) (p_inv_visc_fvf.size ());
        BS_ASSERT (gor_inv_fvf.size () && gor_inv_visc.size () && gor_inv_visc_fvf.size ()) (gor_inv_fvf.size ()) (gor_inv_visc.size ()) (gor_inv_visc_fvf.size ());
        BS_ASSERT (gas_oil_ratio.size () && d_gas_oil_ratio.size () && old_gas_oil_ratio.size ()) (gas_oil_ratio.size ()) (d_gas_oil_ratio.size ()) (old_gas_oil_ratio.size ());
        if (!inv_fvf.size () || !p_inv_fvf.size () || !inv_visc.size () || !p_inv_visc.size () || !inv_visc_fvf.size () || !p_inv_visc_fvf.size ())
          {
            throw bs_exception ("py_pvt_oil::calc", "invalid argument");
          }
        if (!gor_inv_fvf.size () || !gor_inv_visc.size () || !gor_inv_visc_fvf.size ())
          {
            throw bs_exception ("py_pvt_oil::calc", "invalid argument (gor group)");
          }
        if (!gas_oil_ratio.size () || !d_gas_oil_ratio.size () || !old_gas_oil_ratio.size ())
          {
            throw bs_exception ("py_pvt_oil::calc", "invalid argument (gas_oil_ratio group)");
          }

        // TODO: BUG
        index_t i = n * n_phases + 2;
        get_spx (this)->calc (is_g, main_var, p, gor, 
          &inv_fvf [i], &p_inv_fvf[i], &gor_inv_fvf[n], 
          &inv_visc[i], &p_inv_visc[i], &gor_inv_visc[n], 
          &inv_visc_fvf[i], &p_inv_visc_fvf[i], &gor_inv_visc_fvf[n],
          &gas_oil_ratio[n], &d_gas_oil_ratio[n], drsdt, dt, old_gas_oil_ratio[n]);
      }
    };


    void
    py_export_pvt ();

  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
#endif  // #ifndef PY_BS_PVT_H_
