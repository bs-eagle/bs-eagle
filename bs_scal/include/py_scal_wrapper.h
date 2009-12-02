/**
 * \file py_scal_wrapper.h
 * \brief python wrapper for scal_*
 * \author Sergey Miryanov
 * \date 21.05.2008
 * */
#ifndef PY_BS_SCAL_WRAPPER_H_
#define PY_BS_SCAL_WRAPPER_H_

#include "scal_3p.h"
#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"
#include "jfunction.h"

namespace blue_sky
  {
  namespace python
    {

    ////////////////////////////////////////////////////////////////////////////
    //class py_spof_data_source_wrapper : public spof_data_source
    //  {
    //  public:

    //    py_spof_data_source_wrapper (PyObject *obj)
    //        : obj (obj)
    //    {
    //      Py_INCREF (obj);
    //    }
    //    ~py_spof_data_source_wrapper ()
    //    {
    //      Py_DECREF (obj);
    //    }

    //    double Sp () const;
    //    double So () const;
    //    double Krp () const;
    //    double Krop () const;
    //    double Pcp () const;

    //    void first ();
    //    void next ();
    //    bool end ();

    //  private:

    //    PyObject *obj;
    //  };

    //class py_spfn_data_source_wrapper : public spfn_data_source
    //  {
    //  public:

    //    py_spfn_data_source_wrapper  (PyObject *obj)
    //        : obj (obj)
    //    {
    //      Py_INCREF (obj);
    //    }
    //    ~py_spfn_data_source_wrapper ()
    //    {
    //      Py_DECREF (obj);
    //    }

    //    double Sp () const;
    //    double Krp () const;
    //    double Pcp () const;

    //    void first ();
    //    void next ();
    //    bool end ();

    //  private:

    //    PyObject *obj;
    //  };

    //class py_sof3_data_source_wrapper : public sof3_data_source
    //  {
    //  public:

    //    py_sof3_data_source_wrapper (PyObject *obj)
    //        : obj (obj)
    //    {
    //      Py_INCREF (obj);
    //    }

    //    ~py_sof3_data_source_wrapper ()
    //    {
    //      Py_DECREF (obj);
    //    }

    //    double So () const;
    //    double Krow () const;
    //    double Krog () const;

    //    void first ();
    //    void next ();
    //    bool end ();

    //  private:

    //    PyObject *obj;
    //  };
    ////////////////////////////////////////////////////////////////////////////

    template <typename strategy_t>
    class py_jfunction : public py_objbase
      {
      public:

        typedef py_objbase            base_t;
        typedef jfunction<strategy_t> wrapped_t;

        typedef typename strategy_t::item_t item_t;

        py_jfunction ()
            : base_t (wrapped_t::bs_type ())
        {

        }

        void init (JFUNC_PERM_TYPE_ENUM perm_type_)
        {
          get_spx (this)->init (perm_type_);
        }

        item_t get_st_phase () const
          {
            return get_spx (this)->st_phase;
          }
        item_t get_alpha () const
          {
            return get_spx (this)->alpha;
          }
        item_t get_beta () const
          {
            return get_spx (this)->beta;
          }

        void set_st_phase (item_t st_phase)
        {
          this->get_spx (this)->st_phase = st_phase;
        }
        void set_alpha (item_t alpha)
        {
          this->get_spx (this)->alpha = alpha;
        }
        void set_beta (item_t beta)
        {
          this->get_spx (this)->beta = beta;
        }

        bool valid () const
          {
            return get_spx (this)->valid ();
          }

      };

    template <typename strategy_t>
    class py_scale_array_holder : public py_objbase
      {
      public:

        typedef py_objbase                              base_t;
        typedef typename strategy_t::item_t             item_t;
        typedef scale_array_holder <strategy_t>         scale_array_holder_t;
        typedef smart_ptr <scale_array_holder_t, true>  sp_scale_array_holder_t;
        typedef scale_array_holder_t                    wrapped_t;

        py_scale_array_holder ()
            : base_t (wrapped_t::bs_type ())
        {

        }

        py_scale_array_holder (sp_scale_array_holder_t holder)
            : base_t (holder)
        {

        }

        void insert_socr (const shared_vector <item_t> &socr)
        {
          get_spx (this)->insert_socr (socr);
        }
        void insert_scr (const shared_vector <item_t> &scr)
        {
          get_spx (this)->insert_scr (scr);
        }
        void insert_su (const shared_vector <item_t> &su)
        {
          get_spx (this)->insert_su (su);
        }
        void insert_sl (const shared_vector <item_t> &sl)
        {
          get_spx (this)->insert_sl (sl);
        }
      };

    template <typename strategy_t>
    class py_scal_data_holder : public py_objbase
      {
      public:

        typedef py_objbase													base_t;
        typedef scal_2p_data_holder <strategy_t>		scal_2p_data_holder_t;
        typedef scal_2p_data_holder_t								wrapped_t;
        typedef smart_ptr <wrapped_t, true>					sp_scal_data_holder_t;

        py_scal_data_holder ();
        py_scal_data_holder (sp_scal_data_holder_t holder);

#ifdef _DEBUG

        void
        save_raw_data (const std::string &filename);

        void
        save_data (const std::string &filename);

#endif
      };

    //////////////////////////////////////////////////////////////////////////
    //template <typename strategy_t>
    //class py_test_model : public py_objbase
    //  {
    //  public:

    //    typedef py_objbase                       base_t;
    //    typedef scal::test_model <strategy_t>    wrapped_t;

    //    py_test_model ()
    //        : base_t (wrapped_t::bs_type ())
    //    {

    //    }

    //    void init (int n_phases, int cell_count)
    //    {
    //      get_lspx (this)->init (n_phases, cell_count);
    //    }

    //  };
    //////////////////////////////////////////////////////////////////////////

    template <typename strategy_t>
    class py_scal_3p : public py_objbase
      {
      public:

        typedef py_objbase                          base_t;
        typedef typename strategy_t::index_t        index_t;
        typedef typename strategy_t::item_t         item_t;
        typedef typename strategy_t::index_array_t  index_array_t;
        typedef typename strategy_t::item_array_t   item_array_t;
        typedef scal_3p <strategy_t>                scal_3p_t;
        typedef scal_3p_t                           wrapped_t;
        typedef smart_ptr <scal_3p_t, true>         sp_scal_3p_t;

      public:

        py_scal_3p ()
        : base_t (wrapped_t::bs_type ())
        {
        }

        py_scal_3p (const sp_scal_3p_t &scal)
        : base_t (scal)
        {

        }

        void
        get_relative_perm (index_t cell_index, const item_array_t &saturation, const index_array_t &sat_regions, item_array_t &relative_perm, item_array_t &s_deriv_relative_perm) const
        {
          get_spx (this)->get_relative_perm (cell_index, saturation, sat_regions, relative_perm, s_deriv_relative_perm);
        }

        void
        get_capillary (index_t cell_index, const item_array_t &saturation, const index_array_t &sat_regions, const item_array_t &perm, const item_array_t &poro, item_array_t &cap, item_array_t &s_deriv_cap) const
        {
          get_spx (this)->get_capillary (cell_index, saturation, sat_regions, perm, poro, cap, s_deriv_cap);
        }

        py_scale_array_holder<strategy_t> get_water_scale ()
        {
          return py_scale_array_holder<strategy_t> (get_spx (this)->get_water_scale());
        }
        py_scale_array_holder<strategy_t> get_gas_scale ()
        {
          return py_scale_array_holder<strategy_t> (get_spx (this)->get_gas_scale ()); // BUG:
        }

        py_scal_data_holder<strategy_t>
        get_water_data ()
        {
          return py_scal_data_holder<strategy_t> (get_spx (this)->get_water_data ());
        }
        py_scal_data_holder<strategy_t>
        get_gas_data ()
        {
          return py_scal_data_holder<strategy_t> (get_spx (this)->get_gas_data ());
        }

        //void set_water_jfunction (py_jfunction jfunc)
        //{
        //  get_lspx (this)->set_water_jfunction (jfunc.get_spx <jfunction> ());
        //}
        //void set_gas_jfunction (py_jfunction jfunc)
        //{
        //  get_lspx (this)->set_gas_jfunction (jfunc.get_spx <jfunction> ());
        //}

        //void set_n_phases (int n_phases)
        //{
        //  get_lspx (this)->set_n_phases (n_phases);
        //}
        //void set_phases (/*PHASE_ENUM*/int phases)
        //{
        //  get_lspx (this)->set_phases ((PHASE_ENUM)phases);
        //}
        //void set_rpo_model (RPO_MODEL_ENUM rpo_model)
        //{
        //  get_lspx (this)->set_rpo_model (rpo_model);
        //}

        //void set_test_model (py_test_model<strategy_t> tm)
        //{
        //  test_model_ = tm.get_spx <typename py_test_model<strategy_t>::wrapped_t> ();
        //}

        //void set_phase_d (const scal_3p::phase_d_t &phase_d)
        //{
        //  get_lspx (this)->set_phase_d (phase_d);
        //}
        //void set_sat_d (const scal_3p::phase_d_t &sat_d)
        //{
        //  get_lspx (this)->set_sat_d (sat_d);
        //}
      };

    void py_export_scal ();


  } // namespace python
} // namespace blue_sky


#endif // #ifndef PY_BS_SCAL_WRAPPER_H_
