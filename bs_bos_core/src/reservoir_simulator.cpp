/**
 *       \file  reservoir_simulator.cpp
 *      \brief  reservoir_simulator implementation
 *     \author  Nikonov Max
 *       \date  10.11.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "event_base.h"
#include "event_manager.h"
#include "reservoir_simulator.h"
#include "reservoir.h"
#include "data_storage_interface.h"
#include "jacobian.h"
#include "trans_multipliers.h"
#include "main_loop_calc.h"
#include "facility_manager.h"
#include "keyword_manager.h"
#include "strategy_name.h"

namespace blue_sky
{

  /**
   * \brief  'default' ctor for reservoir_simulator
   * \param  param additional params for reservoir_simulator
   * */
  template <class strategy_t>
  reservoir_simulator<strategy_t>::reservoir_simulator (bs_type_ctor_param)
      : bs_refcounter ()
      , bs_node (bs_node::create_node()) //new typename this_t::mstatus_traits ()))
      , dm (give_kernel::Instance().create_object (dm_t::bs_type ()))
      , em (give_kernel::Instance().create_object (em_t::bs_type ()))
      , cm (give_kernel::Instance().create_object (calc_model_t::bs_type ()))
      , mesh (0)
      , reservoir_ (give_kernel::Instance().create_object (reservoir_t::bs_type ()))
      , facility_storage_ (give_kernel::Instance().create_object (facility_storage_t::bs_type ()))
      , jacobian_ (give_kernel::Instance().create_object (jacobian_t::bs_type ()))
      , keyword_manager_ (BS_KERNEL.create_object (keyword_manager_t::bs_type ()))
      , mloop (0)
      , reservoir_simulator_events_init_ (this)
  {
    this->add_signal (BS_SIGNAL_RANGE (reservoir_simulator));
    keyword_manager_->init ();//em);

    //bs_node::insert (bs_link::create (dm, "data_manager"), false);
    //bs_node::insert (bs_link::create (em, "event_manager"), false);
    //bs_node::insert (bs_link::create (cm, "calc_model"), false);

    //bs_node::insert (bs_link::create (reservoir_, "reservoir"), false);
    //bs_node::insert (bs_link::create (jacobian_, "jacobian"), false);
    //bs_node::insert (bs_link::create (facility_storage_, "facility_storage"), false);
  }

  /**
   * \brief  copy-ctor for reservoir_simulator
   * \param  src source copy of reservoir_simulator
   * */
  template <class strategy_t>
  reservoir_simulator<strategy_t>::reservoir_simulator (const this_t &src)
      : bs_refcounter ()
      , bs_node (src)
      , dm (src.dm)
      , em (src.em)
      , cm (src.cm)
      , mesh (src.mesh)
      , reservoir_ (src.reservoir_)
      , facility_storage_ (src.facility_storage_)
      , jacobian_ (src.jacobian_)
      , mloop (0)
      , reservoir_simulator_events_init_ (src.reservoir_simulator_events_init_)
  {
  }

  template <class strategy_t>
  reservoir_simulator<strategy_t>::~reservoir_simulator ()
  {
  }

  template <class strategy_t>
  void reservoir_simulator<strategy_t>::set_mesh (const sp_mesh_iface_t &src)
  {
    this->mesh = src;
  }

  /**
   * \brief  actions that should be executed before read of data
   * \param  em pointer to event_manager instance
   * */
  template <typename sp_em_t>
  void 
  pre_read (const sp_em_t &/*em*/)
  {
  }

  /**
   * \brief  actions that should be executed after read of data
   * \param  em pointer to event_manager instance
   * \details if event_list is empty inserts two events otherwise
   *          if last element of event list contains events we add 
   *          another empty element with date 30 days later than 
   *          previous
   * */
  template <typename em_t>
  void 
  post_read (const smart_ptr <em_t, true> &em)
  {
    //if last element of event list contains events we add another empty element with date 30 days later than previous
    if (em->event_list.size ())
      {
        if (!(*(--em->event_list.end())).second.empty())
          {
            std::list <typename em_t::sp_event_base> tl;
            em->event_list.insert(std::make_pair((*(--em->event_list.end())).first + boost::posix_time::hours(30*24),tl));
          }
      }
    else
      {
        // TODO: BUG:
        BS_ASSERT (false && "em->event_list.empty () == true");

        std::list <typename em_t::sp_event_base> tl;
        em->event_list.insert(std::make_pair(boost::gregorian::from_us_string ("01-01-1970"), tl));
        em->event_list.insert(std::make_pair(boost::gregorian::from_us_string ("01-02-1970"), tl));
      }
  }

  /**
   * \brief  reads model file and process keywords
   * \param  filename name of model file
   * \param  key_handlers pointer to keyword_manager instance
   * \param  params params that will be passed to each keyword_handler
   * */
  template <typename keyword_manager_t>
  void
  read_keyword_file (const std::string &filename, const smart_ptr <keyword_manager_t, true> &key_handlers, typename keyword_manager_t::keyword_params_t &params)
  {
    const typename keyword_manager_t::sp_reader_t &l_reader (params.reader);
    typedef event_manager<typename keyword_manager_t::strategy_type> event_manager_t;
    typedef smart_ptr <event_manager_t, true> sp_event_manager_t;
    sp_event_manager_t em (params.em);
    
    char buf[CHAR_BUF_LEN];
    char key[CHAR_BUF_LEN];
    int flag;
    int len;

    l_reader->init (filename, filename);

    // start of loop for data file reading
    flag = 1;

    for (; flag;)
      {
        // reading keyword
        len = l_reader->read_line (buf, CHAR_BUF_LEN);
        if (len < 0)              // break if EOF
          {
            switch (len)
              {
              case YS_NOTHING_TO_READ:
                break;
              default:
                return;
              }
            //rep->print (LOG_READ_SECTION, LOG_DEBUG, "Finish reading\n");
            BOSOUT (section::read_data, level::low) << "Finish reading" << bs_end;
            break;
          }
        if (sscanf (buf, "%s", key) != 1) // look from buffer keyword
          continue;

        std::string keywrd(key);

        if (key[0] == '/')
          {
            continue;
          }

        if (std::string (key) == "END")
          {
            BOSOUT (section::read_data, level::low) << "Finish reading with END keyword" << bs_end;
            break;
          }

        key_handlers->handle_keyword (keywrd, params);
      }
  }

  template <class strategy_t>
  void reservoir_simulator<strategy_t>::simulate (const std::string &path)
  {
    read_keyword_file_and_init (path);
    main_loop();
  }

  template <typename strategy_t>
  void
  reservoir_simulator <strategy_t>::read_keyword_file_and_init (const std::string &path)
  {
    write_time_to_log init_time ("read model and init reservoir simulator", "");

    keyword_params <strategy_t> params (keyword_manager_, dm->reader, dm->data, em, mesh, cm->ts_params, cm->scal_prop);
    model_filename_ = path;
    pre_read (em); 
    on_pre_read (this);
    read_keyword_file(path, keyword_manager_, params);
    post_read (em);
    on_post_read ();
    //TODO: make keyword handlers able to change pointers on created objects (like mesh)
    mesh = sp_mesh_iface_t (params.mesh, bs_dynamic_cast ());
    if (!mesh)
       {
         bs_throw_exception ("Can't down cast params.mesh");
       }
    init();
  }

  template <class strategy_t>
  const typename reservoir_simulator<strategy_t>::sp_dm_t &
  reservoir_simulator<strategy_t>::get_data_manager () const
  {
    return dm;
  }

  template <class strategy_t>
  const typename reservoir_simulator<strategy_t>::sp_em_t &
  reservoir_simulator<strategy_t>::get_event_manager () const
  {
    return em;
  }

  template <class strategy_t>
  const typename reservoir_simulator<strategy_t>::sp_calc_model_t &
  reservoir_simulator<strategy_t>::get_calc_model () const
  {
    return cm;
  }

  template <class strategy_t>
  const typename reservoir_simulator<strategy_t>::sp_mesh_iface_t &
  reservoir_simulator<strategy_t>::get_mesh() const
  {
    return mesh;
  }

  template <class strategy_t>
  const typename reservoir_simulator<strategy_t>::sp_jacobian_t &
  reservoir_simulator<strategy_t>::get_jacobian () const
  {
    return jacobian_;
  }

  template <typename strategy_t>
  const typename reservoir_simulator <strategy_t>::sp_reservoir_t &
  reservoir_simulator <strategy_t>::get_reservoir () const
  {
    return reservoir_;
  }

  /**
   * \brief  checks saturation saturation
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  template <class strategy_t>
  void 
  check_sat (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &data)
  {
    typedef typename strategy_t::index_t        index_t;
    typedef typename strategy_t::item_t         item_t;
    typedef typename strategy_t::index_array_t  index_array_t;
    typedef typename strategy_t::item_array_t   item_array_t;

    std::ostringstream out_s;

    index_t nb = mesh->get_n_active_elements ();
    const index_array_t &original_element_num = mesh->get_int_to_ext();

    if (!data->init_section)
      {
        if (!(*data->d_map)[PRESSURE].array.size ())
          {
            bs_throw_exception ("check_sat: Initial pressure is not specified");
          }
        if (!(*data->d_map)[SWAT].array.size ())
          {
            bs_throw_exception ("check_sat: Initial saturation is not specified");
          }
        if (!(*data->d_map)[SOIL].array.size ())
          {
            BOSWARN (section::check_data, level::warning) << "SOIL is not initialized" << bs_end;
          }
      }

    if (!data->init_section)
      {
        array_float16_t swat = ((*data->d_map)[SWAT].array);
        array_float16_t soil = ((*data->d_map)[SOIL].array);

        for (index_t i = 0; i < nb; ++i)
          {
            index_t i_m = original_element_num[i];

            if (soil.size () && (soil[i_m] < -1.e-12 || soil[i_m] > 1 + 1.e-5))
              {
                bs_throw_exception (boost::format ("soil[%d] = %f is out of range") % i_m % soil[i_m]);
              }
            if (swat.size () && (swat[i_m] < -1.e-12 || swat[i_m]> 1 + 1.e-5))
              {
                bs_throw_exception (boost::format ("swat[%d] = %f is out of range") % i_m % swat[i_m]);
              }
            if (soil.size () && swat.size () && (soil[i_m] + swat[i_m] > 1 + 1.e-5))
              {
                bs_throw_exception (boost::format ("Sum of water saturation and oil saturation for cell %d = %f is out of range")
                  % i_m % (soil[i_m] + swat[i_m]));
              }
          }
      }
  }

  /**
   * \brief  checks permx, permy, permz
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  template <class strategy_t>
  void 
  check_perm (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &data)
  {
    typedef typename strategy_t::index_t index_t;

    const smart_ptr <idata, true> &ldata (data);
    std::ostringstream out_s;
    index_t nb = mesh->get_n_active_elements ();

    if (!nb)
      {
        bs_throw_exception ("All dimensions should be greate than 0. Number of active elements = 0.");
      }
    if (!ldata->contain ("PERMX") || !ldata->contain ("PERMY") || !ldata->contain ("PERMZ"))
      {
        if (!ldata->contain ("PERMX"))
          {
            BOSERR (section::check_data, level::error) << "PERMX is not specified" << bs_end;
          }
        if (!ldata->contain ("PERMY"))
          {
            BOSERR (section::check_data, level::error) << "PERMY is not specified" << bs_end;
          }
        if (!ldata->contain ("PERMZ"))
          {
            BOSERR (section::check_data, level::error) << "PERMZ is not specified" << bs_end;
          }

        bs_throw_exception ("PERMX, PERMY or PERMZ are not specified");
      }
  }

  /**
   * \brief  checks poro, ntg and multpv
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  template <class strategy_t>
  void 
  check_poro_ntg_multpv (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &data)
  {
    typedef typename strategy_t::index_t index_t;
    typedef typename strategy_t::item_t item_t;
    typedef typename strategy_t::index_array_t index_array_t;
    typedef typename strategy_t::item_array_t item_array_t;

    const index_array_t &original_element_num = mesh->get_int_to_ext();
    index_t nb = mesh->get_n_active_elements ();
    if (!nb)
      {
        bs_throw_exception ("All dimensions should be greate than 0. Number of active elements = 0.");
      }

    array_float16_t poro    = data->get_float_non_empty_array ("PORO");
    array_float16_t ntg     = data->get_float_non_empty_array ("NTG");
    array_float16_t multpv  = data->get_float_array ("MULTPV");

    for (index_t i = 0; i < nb; ++i)
      {
        index_t i_m = original_element_num[i];

        if (poro[i_m] < (item_t) -COMP_EPSILON || poro[i_m] > (item_t) (1 + COMP_EPSILON))
          {
            bs_throw_exception (boost::format ("PORO for node [%d] = %f is out of range")
              % i % poro[i_m]);
          }

        if (poro[i_m] < float16_t (COMP_EPSILON))
          {
            BOSWARN (section::check_data, level::warning) 
              << boost::format ("PORO for node [%d] = %f is too small") % i % poro[i_m] 
              << bs_end;
            poro[i_m] = float16_t (COMP_EPSILON);
          }

        if (ntg[i_m] < float16_t (-COMP_EPSILON))
          {
            bs_throw_exception (boost::format ("NTG for node [%d] = %f is out of range")
              % i % ntg[i_m]);
          }

        if (ntg[i_m] < float16_t (COMP_EPSILON))
          {
            BOSWARN (section::check_data, level::warning) 
              << boost::format ("NTG for node [%d] = %f is too small") % i % ntg[i_m] 
              << bs_end;
            ntg[i_m] = float16_t (COMP_EPSILON);
          }

        if (multpv.size ())
          {
            if (multpv[i_m] < float16_t (-COMP_EPSILON))
              {
                bs_throw_exception (boost::format ("MULTPV for node [%d] = %f is out of range")
                  % i % multpv[i_m]);
              }

            if (multpv[i_m] < float16_t (COMP_EPSILON))
              {
                BOSWARN (section::check_data, level::warning) 
                  << boost::format ("MULTPV for node [%d] = %f is too small") % i % multpv[i_m] 
                  << bs_end;
                multpv[i_m] = float16_t (COMP_EPSILON);
              }
          }

      }
  }

  /**
   * \brief  checks nums - eqlnum, satnum, pvtnum, fipnum
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  template <typename strategy_t>
  void 
  check_num (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &data)
  {
    typedef typename strategy_t::index_t index_t;
    typedef typename strategy_t::index_array_t index_array_t;

    index_t nb = mesh->get_n_active_elements ();

    //Check equil regions number
    if (data->i_map->contain(EQLNUM))
      {
        const array_uint8_t &eqlnum = ((*data->i_map)[EQLNUM].array);

        if (eqlnum.size ())
          {
            for (index_t i = 0; i < nb; ++i)
              {
                index_t i_m = mesh->get_int_to_ext()[i];
                if (eqlnum[i_m] <= 0 || eqlnum[i_m] > data->eql_region)
                  {
                    bs_throw_exception (boost::format ("Equlibrium region number for node [%d] = %d is out of range")
                      % i % eqlnum [i_m]);
                  }
              }
          }
      }

    //Check saturation regions number
    if (data->i_map->contain(SATNUM))
      {
        const array_uint8_t &satnum = tools::get_non_empty (((*data->i_map)[SATNUM].array));

        for (index_t i = 0; i < nb; ++i)
          {
            index_t i_m = mesh->get_int_to_ext()[i];
            if (satnum[i_m] <= 0 || satnum[i_m] > data->sat_region)
              {
                bs_throw_exception (boost::format ("Saturation region number for node [%d] = %d is out of range")
                  % i % satnum [i_m]);
              }
          }
      }

    //Check pvt regions number
    if (data->i_map->contain(PVTNUM))
      {
        const array_uint8_t &pvtnum = tools::get_non_empty (((*data->i_map)[PVTNUM].array));

        for (index_t i = 0; i < nb; ++i)
          {
            index_t i_m = mesh->get_int_to_ext()[i];
            if (pvtnum[i_m] <= 0 || pvtnum[i_m] > data->pvt_region)
              {
                bs_throw_exception (boost::format ("PVT region number for node [%d] = %d is out of range")
                  % i % pvtnum [i_m]);
              }
          }
      }

    //Check fip regions number
    if (data->i_map->contain(FIPNUM))
      {
        const array_uint8_t &fipnum = tools::get_non_empty ((*data->i_map)[FIPNUM].array);

        for (index_t i = 0; i < nb; ++i)
          {
            index_t i_m = mesh->get_int_to_ext()[i];
            if (fipnum[i_m] <= 0 || fipnum[i_m] > data->fip_region)
              {
                bs_throw_exception (boost::format ("FIP region number for node [%d] = %d is out of range")
                  % i % fipnum [i_m]);
              }
          }
      }

  }

  /**
   * \brief  checks rock and p_ref
   * \param  data pointer to data storage (pool)
   * \return 
   * */
  template <class strategy_t>
  void 
  check_rock (const smart_ptr <idata, true> &data)
  {
    typedef typename strategy_t::index_t index_t;
    std::ostringstream out_s;
    if (!data->rock.size ())
      {
        bs_throw_exception ("ROCK is not specified");
      }

    for (size_t i = 0; i < data->pvt_region; ++i)
      {
        if (data->rock[i] < 0)
          {
            bs_throw_exception (boost::format ("Rock compressibility in region %d = %f is out of range")
              % (i + 1) % data->rock [i]);
          }
        if (data->p_ref[i] <= 1.e-14)
          {
            bs_throw_exception (boost::format ("Reference pressure for rock compressibility in region %d = %f is out of range")
              % (i + 1) % data->p_ref[i]);
          }
        // Both compressibilities (rock and water) cannot be 0 both
        if (data->rock[i] < COMP_EPSILON && data->pvtw[i].main_data_[2] < COMP_EPSILON)   // 2 = DATA_GPR in pvt_water
          {
            // TODO: BUG: data->pvtw[i].main_data_[2]
            bs_throw_exception (boost::format ("In region %d rock compressibility = %f and water compressibility = %f are both about 0")
              % (i + 1) % data->rock[i] % data->pvtw[i].main_data_[2]);
          }
      }
  }


  /**
   * \brief  checks geometry (mesh), redirect call to mesh
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  template <class strategy_t>
  void 
  check_geometry (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &data)
  {
    typedef typename strategy_t::index_t index_t;
    typedef typename strategy_t::item_t item_t;
    typedef typename strategy_t::index_array_t index_array_t;
    typedef typename strategy_t::item_array_t item_array_t;

    std::ostringstream out_s;

    index_t nb = data->dimens.nx * data->dimens.ny * data->dimens.nz;
    
    if (!nb)
      {
        bs_throw_exception ("All dimensions should be greater than 0. Number of active elements = 0.");
      }
    
    
    
    /*
    // In case of restart we got all geometry from mesh restored from
    // the restart file. So we skip this test.
    if (data->restart)
      return;
    */
    
    mesh->check_data();
  }

  /**
   * \brief  checks equil
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  template <class strategy_t>
  void 
  check_equil (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &data)
  {
    typedef typename strategy_t::index_t index_t;
    typedef typename strategy_t::item_t item_t;
    typedef typename strategy_t::index_array_t index_array_t;

    std::ostringstream out_s;
    index_t N_eq = 2; // number of regions (eql + pvt = 2)
    index_t nb = mesh->get_n_active_elements ();

    index_t eq_reg = N_eq * data->eql_region;
    data->equil_regions.resize (eq_reg, 0);

    array_uint8_t eqlnum;
    array_uint8_t pvtnum;
    if (!data->i_map->contain(EQLNUM))
      {
        data->i_map->create_item (EQLNUM, &i_pool_sizes[ARRAY_POOL_TOTAL * EQLNUM], i_pool_default_values[EQLNUM]);
        eqlnum = tools::get_non_empty ((*data->i_map)[EQLNUM].array);
        eqlnum.assign (i_pool_default_values[EQLNUM] + 1);
      }
    if (!data->i_map->contain(PVTNUM))
      {
        data->i_map->create_item(PVTNUM, &i_pool_sizes[ARRAY_POOL_TOTAL * PVTNUM], i_pool_default_values[PVTNUM]);
        pvtnum = tools::get_non_empty ((*data->i_map)[PVTNUM].array);
        pvtnum.assign (i_pool_default_values[PVTNUM] + 1);
      }

    eqlnum = tools::get_non_empty ((*data->i_map)[EQLNUM].array);
    pvtnum = tools::get_non_empty ((*data->i_map)[PVTNUM].array);

    //check correspondence of pvt-regions to eql-regions:
    //one eql-region must contain only one pvt-region
    for (index_t i = 0; i < nb; ++i)
      {
        index_t i_m = mesh->get_int_to_ext ()[i];
        eq_reg = (eqlnum[i_m] - 1) * N_eq;
        if (data->equil_regions[eq_reg] == 0)
          {
            // set values of regions for current equil region
            data->equil_regions[eq_reg] = 1;
            data->equil_regions[eq_reg + 1] = pvtnum[i_m] - 1;   // set pvt region for equil region
          }
        else
          {
            if (pvtnum[i_m] != data->equil_regions[eq_reg + 1] + 1)
              {
                bs_throw_exception (boost::format ("Incorrect PVT region %d for equilibrium region %d in grid block %d")
                  % pvtnum [i_m] % eqlnum [i_m] % i_m);
              }
          }
      }
  }

  /**
   * \brief  checks volume
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  template <class strategy_t>
  void 
  check_volume (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &/*data*/)
  {
    typedef typename strategy_t::index_t index_t;
    typedef typename strategy_t::item_array_t item_array_t;
    const smart_ptr <rs_smesh_iface <strategy_t>, true> s_mesh(mesh, bs_dynamic_cast ());
    std::ostringstream out_s;
    index_t x, y, z;
    
    const item_array_t &mesh_volumes = mesh->get_volumes ();
    for (index_t i = 0, cnt = mesh->get_n_active_elements(); i < cnt; ++i)
      {
        if (mesh_volumes[i] < COMP_EPSILON)
          {
            //TODO: wrap mesh element identification to single function for all mesh types
            s_mesh->get_element_int_to_ijk (i, x, y, z);

            bs_throw_exception (boost::format ("Volume for grid block [%d: %d, %d, %d] = %f is too small")
              % i % x % y % z % mesh_volumes[i]);
          }
      }
  }


  /**
   * \todo write properly comment
   * */
  template <typename T, typename sp_idata_t, typename sp_mesh_iface_t>
  bool
  check5 (T function, const sp_idata_t &ldata, const sp_mesh_iface_t &mesh, const char *description, bool check_init)
  {
    try
      {
        BOSOUT (section::check_data, level::medium) << description << bs_end;
        if (check_init && !ldata->init_section)
          {
            BOSOUT (section::check_data, level::medium) << " ...[ NOT INITED ]" << bs_end;
          }
        else
          {
            function (mesh, ldata);
            return true;
          }
      }
    catch (bs_exception &e)
      {
        BOSOUT (section::check_data, level::medium) << " ...[FAIL]" << bs_end;
        BOSERR (section::check_data, level::error) << e.what () << bs_end;
      }

    return false;
  }

  /**
   * \todo write properly comment
   * */
  template <typename T, typename sp_idata_t>
  bool
  check3 (T function, const sp_idata_t &ldata, const char *description)
  {
    try
      {
        BOSOUT (section::check_data, level::medium) << description << bs_end;

        function (ldata);
        return true;
      }
    catch (bs_exception &e)
      {
        BOSOUT (section::check_data, level::medium) << " ...[FAIL]" << bs_end;
        BOSERR (section::check_data, level::error) << e.what () << bs_end;
      }

    return false;
  }

  /**
   * \brief  checks all data
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * */
  template <typename strategy_t>
  void
  check_data (const smart_ptr <rs_mesh_iface <strategy_t>, true> &mesh, const smart_ptr <idata, true> &data)
  {
    std::ostringstream out_s;

    BOSOUT (section::check_data, level::medium)
    << "\n***************************************************************\n"
    << "                            CHECK INPUT DATA\n"
    << "***************************************************************" << bs_end;

    size_t count = 0;
    count += check5 (check_num <strategy_t>,              data, mesh,  "Region numbers check... ", true);
    count += check5 (check_sat <strategy_t>,              data, mesh,  "Saturation check... ", false);
    count += check5 (check_geometry <strategy_t>,         data, mesh,  "Geometry check... ", false);
    count += check5 (check_perm <strategy_t>,             data, mesh,  "Permability check... ", false);
    count += check5 (check_poro_ntg_multpv <strategy_t>,  data, mesh,  "Porosity and NTG, MULTPV check check... ", false);
    count += check5 (check_equil <strategy_t>,            data, mesh,  "Equil region numbers check... ", true);
    count += check5 (check_volume <strategy_t>,           data, mesh,  "Volume check... ", false);
    count += check3 (check_rock <strategy_t>,             data,        "Rock check... ");

    //BOSOUT (section::check_data, level::medium) << "Well data check... ";
    //if (!(r = check_wfrictn_keyword_connection_data (msh)))
    //BOSOUT (section::check_data, level::medium) << " ...[ OK ]" << bs_end;
    //else
    //BOSOUT (section::check_data, level::medium) << " ...[FAIL] - error: " << r << bs_end;

    //if (count < 8)
    //  throw bs_exception ("data_manager", "one of check failed");

    // TODO:
    //BOSOUT << output_time;
  }

  template <class strategy_t>
  void reservoir_simulator<strategy_t>::init ()
  {
    std::ostringstream out_s;

    const typename calc_model_t::sp_jacobian_matrix_t &jmatrix (jacobian_->get_jmatrix ());

    dm->check_arrays_for_inactive_blocks();

    mesh->init_props (dm->data);
    mesh->set_darcy (cm->internal_constants.darcy_constant);

    //!TODO:output init_ext_to_int(splicing)
    mesh->init_ext_to_int();

    if (mesh->get_n_active_elements() <= 0)
      {
        bs_throw_exception ("Error: No active cells!");
      }

    //////////////////////////////////////
    // Check input data before units conversion for correct reporting
    // of user errors
    check_data (mesh, dm->data);

    cm->init_main_arrays(dm->data, mesh);
    cm->init_calcul_arrays (dm->data, mesh);
    cm->init_jacobian (jacobian_, mesh);

    // calculate planes geometric transmissibility
    cm->rock_grid_prop->init_planes_trans (mesh->get_n_active_elements (), mesh->get_volumes (), cm->ts_params, cm->internal_constants);

    typename mesh_iface_t::index_array_t boundary_array;
    std::string flux_conn_name = "bs_flux_connections_" + tools::strategy_name <strategy_t>::name ();
    typename mesh_iface_t::sp_flux_conn_iface_t flux_conn (BS_KERNEL.create_object(flux_conn_name), bs_dynamic_cast ());
    if (!flux_conn)
      {
       bs_throw_exception (boost::format ("Can`t create %s!") % flux_conn_name); 
      }
    
    index_t N_block_size = cm->n_phases;
    index_t N_blocks = mesh->get_n_active_elements ();

    mesh->build_jacobian_and_flux_connections (jmatrix->get_regular_matrix(), flux_conn, boundary_array);
    jmatrix->get_regular_matrix()->init_values(N_block_size);
    assign (jmatrix->get_regular_matrix ()->get_values(), N_block_size * N_block_size * (2* mesh->get_n_connections() + mesh->get_n_active_elements()), item_t (0));
    jmatrix->m_array = flux_conn->get_matrix_block_idx_minus();
    jmatrix->p_array = flux_conn->get_matrix_block_idx_plus ();
    jmatrix->trns_matrix = flux_conn->get_conn_trans();

    jmatrix->get_irregular_matrix ()->init (N_blocks, N_blocks, N_block_size, 0);
  }

  template <typename strategy_t>
  typename reservoir_simulator <strategy_t>::main_loop_calc_t *
  reservoir_simulator <strategy_t>::get_main_loop ()
  {
    if (cm->n_phases == 3)
      return new main_loop_calc <strategy_t, true, true, true> (this);
    else if (cm->is_water () && cm->is_oil ())
      return new main_loop_calc <strategy_t, true, false, true> (this);
    else if (cm->is_gas () && cm->is_oil ())
      return new main_loop_calc <strategy_t, false, true, true> (this);
    else if (cm->is_water ())
      return new main_loop_calc <strategy_t, true, false, false> (this);
    else if (cm->is_gas ())
      return new main_loop_calc <strategy_t, false, true, false> (this);
    else if (cm->is_oil ())
      return new main_loop_calc <strategy_t, false, false, true> (this);
    else
      bs_throw_exception ("Unknown phase model");
  }

  template <class strategy_t>
  void
  reservoir_simulator<strategy_t>::main_loop ()
  {
    on_begin (clock ());
    write_time_to_log lllll("main loop time","tss");

    mloop = get_main_loop ();
    if (!mloop)
      {
        bs_throw_exception ("Can't create main_loop_calc");
      }

    try 
      {
        mloop->ready ();
        mloop->go ();
        mloop->end ();
        mloop = 0;

        on_end (clock ());
      }
    catch (...)
      {
        mloop = 0;

        throw;
      }
  }

  template <typename strategy_t>
  void
  reservoir_simulator <strategy_t>::pre_large_step (const sp_event_base_list_t &event_list)
  {
    mloop->apply_events (event_list);
    reservoir_->pre_large_step (cm, mesh);
  }

  template <typename strategy_t>
  std::string
  reservoir_simulator <strategy_t>::model_filename () const
  {
    return model_filename_;
  }

  // create object
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (reservoir_simulator, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (reservoir_simulator, (class))

  // array map implementation
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (reservoir_simulator <base_strategy_fi>), 1, (bs_node), "reservoir_simulator_fi", "Reservoir simulator fi", "Reservoir simulator float", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (reservoir_simulator <base_strategy_di>), 1, (bs_node), "reservoir_simulator_di", "Reservoir simulator di", "Reservoir simulator double", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (reservoir_simulator <base_strategy_mixi>), 1, (bs_node), "reservoir_simulator_mixi", "Reservoir simulator mixi", "Reservoir simulator mix", false);



  bool
  reservoir_simulator_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, reservoir_simulator<base_strategy_fi>::bs_type ());    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, reservoir_simulator<base_strategy_di>::bs_type ());    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, reservoir_simulator<base_strategy_mixi>::bs_type ());  BS_ASSERT (res);

    return res;
  }
}
