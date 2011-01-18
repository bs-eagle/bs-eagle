#ifndef PY_DATA_MANAGER_H
#define PY_DATA_MANAGER_H

namespace blue_sky
  {

  template <typename strategy_t>
  class data_manager;

  namespace python
  {
    template <typename strategy_t>
    class BS_API_PLUGIN py_idata;

    template <class strategy_t>
    class BS_API_PLUGIN py_data_manager : public py_objbase
      {
      public:
        typedef data_manager<strategy_t>			wrapped_t;
        typedef smart_ptr<wrapped_t, true>		sp_dm_t;
        typedef py_idata <strategy_t>         py_idata_t;

      public:
        //py_data_manager(const py_event_manager_t &);
        py_data_manager(const sp_dm_t &);

        virtual ~py_data_manager();

        py_idata_t get_idata () const;
      };

    void
    py_export_data_manager ();

  }	// namespace python
}	// namespace blue_sky

#endif // PY_DATA_MANAGER_H
