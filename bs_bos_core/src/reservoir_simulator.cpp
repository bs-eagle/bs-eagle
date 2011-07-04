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

namespace blue_sky
{

  /**
   * \brief  'default' ctor for reservoir_simulator
   * \param  param additional params for reservoir_simulator
   * */
  reservoir_simulator::reservoir_simulator (bs_type_ctor_param)
      : bs_refcounter ()
      , bs_node (bs_node::create_node()) //new typename this_t::mstatus_traits ()))
      , hdm_ (give_kernel::Instance().create_object (hdm::bs_type ()))
      , em (give_kernel::Instance().create_object (event_manager::bs_type ()))
      , cm (give_kernel::Instance().create_object (calc_model::bs_type ()))
      , reservoir_ (give_kernel::Instance().create_object (reservoir::bs_type ()))
      , facility_storage_ (give_kernel::Instance().create_object (facility_storage_t::bs_type ()))
      , mloop (0)
      , reservoir_simulator_events_init_ (this)
  {
    this->add_signal (BS_SIGNAL_RANGE (reservoir_simulator));
    hdm_->init ("hdm_model");
    
    //bs_node::insert (bs_link::create (hdm_, "hdm"), false);
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
  reservoir_simulator::reservoir_simulator (const reservoir_simulator &src)
      : bs_refcounter ()
      , bs_node (src)
      , hdm_ (src.hdm_)
      , em (src.em)
      , cm (src.cm)
      , reservoir_ (src.reservoir_)
      , facility_storage_ (src.facility_storage_)
      , jacobian_ (src.jacobian_)
      , mloop (0)
      , reservoir_simulator_events_init_ (src.reservoir_simulator_events_init_)
  {
  }

  reservoir_simulator::~reservoir_simulator ()
  {
  }

  void reservoir_simulator::set_mesh (const sp_mesh_iface_t &mesh)
  {
    // FIXME: is set mesh is good idea?
    this->hdm_->set_mesh (mesh);
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
  void
  read_keyword_file (const std::string &filename, 
    const smart_ptr <hdm_iface, true> &hdm_, 
    const smart_ptr <event_manager, true> &em)
  {
    hdm_iface::sp_reader_t reader = hdm_->get_reader ();
    hdm_iface::sp_km_t keywords = hdm_->get_keyword_manager ();
    keyword_params params (hdm_);
    
    char buf[CHAR_BUF_LEN];
    char key[CHAR_BUF_LEN];
    int flag;
    int len;
    
    write_time_to_log init_time ("Read model", "");

    reader->init (filename, filename);

    // start of loop for data file reading
    flag = 1;
    for (; flag;)
      {
        // reading keyword
        len = reader->read_line (buf, CHAR_BUF_LEN);
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

        keywords->handle_keyword (keywrd, params);
      }
  }

  void reservoir_simulator::simulate (const std::string &path)
  {
    hdm_->init (path);
    read_keyword_file_and_init (path);
    main_loop();
  }

  void
  reservoir_simulator::read_keyword_file_and_init (const std::string &path)
  {
    write_time_to_log init_time ("read model and init reservoir simulator", "");

    //keyword_params params (keyword_manager_, hdm->reader, hdm->data, em, mesh, cm->ts_params, cm->scal_prop);
    model_filename_ = path;
    pre_read (em); 
    on_pre_read (this);
    read_keyword_file(path, hdm_, em);
    post_read (em);
    on_post_read ();
    init();
  }

  reservoir_simulator::sp_hdm_t 
  reservoir_simulator::get_hdm () const
  {
    return hdm_;
  }

  reservoir_simulator::sp_em_t 
  reservoir_simulator::get_event_manager () const
  {
    return em;
  }

  reservoir_simulator::sp_calc_model_t 
  reservoir_simulator::get_calc_model () const
  {
    return cm;
  }

  reservoir_simulator::sp_mesh_iface_t
  reservoir_simulator::get_mesh() const
  {
    return hdm_->get_mesh ();
  }

  reservoir_simulator::sp_jacobian_t 
  reservoir_simulator::get_jacobian () const
  {
    return jacobian_;
  }

  reservoir_simulator::sp_reservoir_t 
  reservoir_simulator::get_reservoir () const
  {
    return reservoir_;
  }

  // FIXME: WTF facility_storage_t
  smart_ptr <reservoir_simulator::facility_storage_t> 
  reservoir_simulator::get_facility_storage () const
  {
    return facility_storage_;
  }

  /**
   * \brief  checks saturation saturation
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  void 
  check_sat (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &data)
  {
    t_long nb = mesh->get_n_active_elements ();
    const spv_long &original_element_num_ = mesh->get_int_to_ext();
    const t_long *original_element_num = &(*original_element_num_)[0];

    if (!data->props->get_i ("init_section"))
      {
        if (!data->contains_fp_array ("PRESSURE"))
          {
            bs_throw_exception ("check_sat: Initial pressure is not specified");
          }
        if (!data->contains_fp_array ("SWAT"))
          {
            bs_throw_exception ("check_sat: Initial saturation is not specified");
          }
        if (!data->contains_fp_array ("SOIL"))
          {
            BOSWARN (section::check_data, level::warning) << "SOIL is not initialized" << bs_end;
          }

        spv_float swat_ = data->get_fp_array ("SWAT");
        spv_float soil_ = data->get_fp_array ("SOIL", false);

        t_float *swat = &(*swat_)[0]; 
        t_float *soil = soil_ && soil_->size () ? &(*soil_)[0] : 0;


        for (t_long i = 0; i < nb; ++i)
          {
            t_long i_m = original_element_num[i];

            if (soil && (soil[i_m] < -1.e-12 || soil[i_m] > 1 + 1.e-5))
              {
                bs_throw_exception (boost::format ("soil[%d] = %f is out of range") % i_m % soil[i_m]);
              }
            if (swat[i_m] < -1.e-12 || swat[i_m]> 1 + 1.e-5)
              {
                bs_throw_exception (boost::format ("swat[%d] = %f is out of range") % i_m % swat[i_m]);
              }
            if (soil && soil[i_m] + swat[i_m] > 1 + 1.e-5)
              {
                bs_throw_exception (boost::format ("Sum of water saturation and oil saturation for cell %d = %f is out of range")
                  % i_m % (soil[i_m] + swat[i_m]));
              }
          }
      }
  }


  /**
   * \brief  checks saturation saturation
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  void 
  check_scal (const smart_ptr <idata, true> &data)
  {
    // TODO: check code 
    typedef idata::scal_vector    scal_vector;
    typedef idata::scal_info      scal_info;

    bool is_water = data->props->get_b ("water_phase");
    bool is_gas   = data->props->get_b ("gas_phase");
    bool is_oil   = data->props->get_b ("oil_phase");
    
    if (is_water && is_gas && is_oil)
      {
        if (!data->swof.front ().main_data_->empty () && !data->sgof.front ().main_data_->empty ())
          {
            data->props->set_i ("scal_family", 0);
          }
        else if (!data->swfn.front ().main_data_->empty () && 
                 !data->sgfn.front ().main_data_->empty () &&
                 !data->sof3.front ().main_data_->empty ())
          {
            data->props->set_i ("scal_family", 1);
          }         
      } 
    else if (is_water && is_oil) 
      {
        if (!data->swof.front ().main_data_->empty ())
          {
            data->props->set_i ("scal_family", 0);
          }
        else if (!data->swfn.front ().main_data_->empty () && 
                 !data->sof2.front ().main_data_->empty ())         
          {
            data->props->set_i ("scal_family", 1);
          }       
      }  
    else if (is_oil && is_gas)
      {
        if (!data->sgof.front ().main_data_->empty ())
          {
            data->props->set_i ("scal_family", 0);
          }
        else if (!data->sgfn.front ().main_data_->empty () && 
                 !data->sof2.front ().main_data_->empty ())         
          {
            data->props->set_i ("scal_family", 1);
          }       
      }  
  }

  /**
   * \brief  checks permx, permy, permz
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  void 
  check_perm (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &data)
  {
    t_long nb = mesh->get_n_active_elements ();
    if (!nb)
      {
        bs_throw_exception ("All dimensions should be greate than 0. Number of active elements = 0.");
      }
    if (!data->contains_fp_array ("PERMX") || !data->contains_fp_array ("PERMY") || !data->contains_fp_array ("PERMZ"))
      {
        if (!data->contains_fp_array ("PERMX"))
          {
            BOSERR (section::check_data, level::error) << "PERMX is not specified" << bs_end;
          }
        if (!data->contains_fp_array ("PERMY"))
          {
            BOSERR (section::check_data, level::error) << "PERMY is not specified" << bs_end;
          }
        if (!data->contains_fp_array ("PERMZ"))
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
  void 
  check_poro_ntg_multpv (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &data)
  {
    t_long nb = mesh->get_n_active_elements ();
    if (!nb)
      {
        bs_throw_exception ("All dimensions should be greate than 0. Number of active elements = 0.");
      }

    spv_double poro_ = data->get_fp_array ("PORO");
    spv_double ntg_ = data->get_fp_array ("NTG");
    spv_double multpv_ = data->get_fp_array ("MULTPV");
    spv_long original_element_num_ = mesh->get_int_to_ext();

    t_long const * original_element_num = &(*original_element_num_)[0];
    t_double * poro = &(*poro_)[0];
    t_double * ntg = &(*ntg_)[0];
    t_double * multpv = multpv_ ? &(*multpv_)[0] : 0;

    for (t_long i = 0; i < nb; ++i)
      {
        t_long i_m = original_element_num[i];

        if (poro[i_m] < (t_double) -COMP_EPSILON || poro[i_m] > (t_double) (1 + COMP_EPSILON))
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

  /**
   * \brief  checks nums - eqlnum, satnum, pvtnum, fipnum
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  void 
  check_num (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &data)
  {
    t_long nb = mesh->get_n_active_elements ();
    spv_long index_map_ = mesh->get_int_to_ext();
    t_long const *index_map = &(*index_map_)[0];

    //Check equil regions number
    if (data->contains_i_array ("EQLNUM"))
      {
        spv_int eqlnum_ = data->get_i_array ("EQLNUM");
        t_int const *eqlnum = &(*eqlnum_)[0];
        t_long eql_region = data->props->get_i ("eql_region");

        // FIXME: check nb and array size
        for (t_long i = 0; i < nb; ++i)
          {
            t_long i_m = index_map[i];
            if (eqlnum[i_m] <= 0 || eqlnum[i_m] > eql_region)
              {
                bs_throw_exception (boost::format ("Equlibrium region number for node [%d] = %d is out of range")
                  % i % eqlnum [i_m]);
              }
          }
      }

    //Check saturation regions number
    if (data->contains_i_array ("SATNUM"))
      {
        spv_int satnum_ = data->get_i_array ("SATNUM");
        t_int const *satnum = &(*satnum_)[0];
        t_long sat_region = data->props->get_i ("sat_region");

        // FIXME: check nb and array size
        for (t_long i = 0; i < nb; ++i)
          {
            t_long i_m = index_map[i];
            if (satnum[i_m] <= 0 || satnum[i_m] > sat_region)
              {
                bs_throw_exception (boost::format ("Saturation region number for node [%d] = %d is out of range")
                  % i % satnum [i_m]);
              }
          }
      }

    //Check pvt regions number
    if (data->contains_i_array ("PVTNUM"))
      {
        spv_int pvtnum_ = data->get_i_array ("PVTNUM");
        t_int const *pvtnum = &(*pvtnum_)[0];
        t_long pvt_region = data->props->get_i ("pvt_region");

        // FIXME: check nb and array size
        for (t_long i = 0; i < nb; ++i)
          {
            t_long i_m = index_map[i];
            if (pvtnum[i_m] <= 0 || pvtnum[i_m] > pvt_region)
              {
                bs_throw_exception (boost::format ("PVT region number for node [%d] = %d is out of range")
                  % i % pvtnum [i_m]);
              }
          }
      }

    //Check fip regions number
    if (data->contains_i_array ("FIPNUM"))
      {
        spv_int fipnum_ = data->get_i_array ("FIPNUM");
        t_int const *fipnum = &(*fipnum_)[0];
        t_long fip_region = data->props->get_i ("fip_region");

        // FIXME: check nb and array size
        for (t_long i = 0; i < nb; ++i)
          {
            t_long i_m = index_map[i];
            if (fipnum[i_m] <= 0 || fipnum[i_m] > fip_region)
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
  void 
  check_rock (const smart_ptr <idata, true> &data)
  {
    if (!data->rock->size ())
      {
        bs_throw_exception ("ROCK is not specified");
      }

    t_float const *rock = &(*data->rock)[0];
    t_float const *p_ref = &(*data->p_ref)[0];
    t_long pvt_region = data->props->get_i ("pvt_region");
    for (t_long i = 0; i < pvt_region; ++i)
      {
        if (rock[i] < 0)
          {
            bs_throw_exception (boost::format ("Rock compressibility in region %d = %f is out of range")
              % (i + 1) % rock [i]);
          }
        if (p_ref[i] <= 1.e-14)
          {
            bs_throw_exception (boost::format ("Reference pressure for rock compressibility in region %d = %f is out of range")
              % (i + 1) % p_ref[i]);
          }
        // Both compressibilities (rock and water) cannot be 0 both
        if (rock[i] < COMP_EPSILON && (*data->pvtw[i].main_data_)[2] < COMP_EPSILON)   // 2 = DATA_GPR in pvt_water
          {
            // TODO: BUG: data->pvtw[i].main_data_[2]
            bs_throw_exception (boost::format ("In region %d rock compressibility = %f and water compressibility = %f are both about 0")
              % (i + 1) % rock[i] % (*data->pvtw[i].main_data_)[2]);
          }
      }
  }


  /**
   * \brief  checks geometry (mesh), redirect call to mesh
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  void 
  check_geometry (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &data)
  {
    if (!mesh->get_n_elements ())
      {
        bs_throw_exception ("Number of elements in mesh = 0");
      }
    
    if (!mesh->get_n_active_elements ())
      {
        bs_throw_exception ("Number of active elements in mesh = 0.");
      }
    
    mesh->check_data();
  }

  /**
   * \brief  checks equil
   * \param  mesh pointer to mesh instance
   * \param  data pointer to data storage (pool)
   * \return may throw exception
   * */
  void 
  check_equil (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &data)
  {
    t_long N_eq = 2; // number of regions (eql + pvt = 2)
    t_long eq_reg = N_eq * data->props->get_i ("eql_region");
    data->equil_regions->init (eq_reg, 0);

    if (!data->contains_i_array ("EQLNUM"))
      {
        data->create_i_array ("EQLNUM", &i_pool_sizes[ARRAY_POOL_TOTAL * EQLNUM], i_pool_default_values[EQLNUM]);
      }
    if (!data->contains_i_array ("PVTNUM"))
      {
        data->create_i_array ("PVTNUM", &i_pool_sizes[ARRAY_POOL_TOTAL * PVTNUM], i_pool_default_values[PVTNUM]);
      }

    spv_int eqlnum_ = data->get_i_array ("EQLNUM");
    spv_int pvtnum_ = data->get_i_array ("PVTNUM");
    t_int const *eqlnum = &(*eqlnum_)[0];
    t_int const *pvtnum = &(*pvtnum_)[0];
    
    t_int *equil_regions = &(*data->equil_regions)[0];

    //check correspondence of pvt-regions to eql-regions:
    //one eql-region must contain only one pvt-region
    spv_long index_map_ = mesh->get_int_to_ext ();
    t_long const *index_map = &(*index_map_)[0];
    t_long nb = mesh->get_n_active_elements ();
    for (t_long i = 0; i < nb; ++i)
      {
        t_long i_m = index_map[i];
        eq_reg = (eqlnum[i_m] - 1) * N_eq;
        if (equil_regions[eq_reg] == 0)
          {
            // set values of regions for current equil region
            equil_regions[eq_reg] = 1;
            equil_regions[eq_reg + 1] = pvtnum[i_m] - 1;   // set pvt region for equil region
          }
        else
          {
            if (pvtnum[i_m] != equil_regions[eq_reg + 1] + 1)
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
  void 
  check_volume (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &/*data*/)
  {
    const smart_ptr <rs_smesh_iface, true> s_mesh(mesh, bs_dynamic_cast ());
    
    const spv_double &mesh_volumes_ = mesh->get_volumes ();
    t_double const *mesh_volumes = &(*mesh_volumes_)[0];
    for (t_long i = 0, cnt = mesh->get_n_active_elements(); i < cnt; ++i)
      {
        if (mesh_volumes[i] < COMP_EPSILON)
          {
            //TODO: wrap mesh element identification to single function for all mesh types
            t_long x, y, z;
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
        if (check_init && !ldata->props->get_i ("init_section"))
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
  void
  check_data (const smart_ptr <rs_mesh_iface, true> &mesh, const smart_ptr <idata, true> &data)
  {
    std::ostringstream out_s;
    
    BOSOUT (section::check_data, level::medium)
    << "\n***************************************************************\n"
    << "                            CHECK INPUT DATA\n"
    << "***************************************************************" << bs_end;

    size_t count = 0;
    count += check5 (check_num,              data, mesh,  "Region numbers check... ", true);
    count += check5 (check_sat,              data, mesh,  "Saturation check... ", false);
    count += check5 (check_geometry,         data, mesh,  "Geometry check... ", false);
    count += check5 (check_perm,             data, mesh,  "Permability check... ", false);
    count += check5 (check_poro_ntg_multpv,  data, mesh,  "Porosity and NTG, MULTPV check ... ", false);
    count += check5 (check_equil,            data, mesh,  "Equil region numbers check... ", true);
    count += check5 (check_volume,           data, mesh,  "Volume check... ", false);
    count += check3 (check_rock,             data,        "Rock check... ");
    count += check3 (check_scal,             data,        "SCAL check... ");
    
    //BOSOUT (section::check_data, level::medium) << "Well data check... ";
    //if (!(r = check_wfrictn_keyword_connection_data (msh)))
    //BOSOUT (section::check_data, level::medium) << " ...[ OK ]" << bs_end;
    //else
    //BOSOUT (section::check_data, level::medium) << " ...[FAIL] - error: " << r << bs_end;

    //if (count < 8)
    //  throw bs_exception ("hdm", "one of check failed");

    // TODO:
    //BOSOUT << output_time;
  }

  void 
  reservoir_simulator::init ()
  {
    hdm_->check_arrays_for_inactive_blocks();

    hdm_->get_mesh ()->init_props (hdm_);
    hdm_->get_mesh ()->set_darcy (cm->internal_constants.darcy_constant);

    //!TODO:output init_ext_to_int(splicing)
    hdm_->get_mesh ()->init_ext_to_int();
    t_long n_active_elements = hdm_->get_mesh ()->get_n_active_elements ();
    if (n_active_elements <= 0)
      {
        bs_throw_exception ("Error: No active cells!");
      }

    // Check input data before units conversion for correct reporting of user errors
    check_data (hdm_->get_mesh (), hdm_->data);

    cm->init_main_arrays(hdm_->get_init_model (), hdm_->get_scal (), hdm_->get_data (), hdm_->get_mesh ());
    cm->init_calcul_arrays (hdm_->data, hdm_->get_mesh ());
    // calculate planes geometric transmissibility
    cm->rock_grid_prop->init_planes_trans (n_active_elements, hdm_->mesh->get_volumes (), cm->ts_params, cm->internal_constants);

    // now jacobian contains flux_connections and boundary array
    jacobian_ = BS_KERNEL.create_object (jacobian::bs_type ());
    jacobian_->init (n_active_elements, cm->n_phases, cm->n_sec_vars);

    hdm_->get_mesh ()->build_jacobian_and_flux_connections (jacobian_->get_matrix ("flux"), 
      jacobian_->get_flux_connections (),
      jacobian_->get_boundary ());
  }

  main_loop_calc_base *
  reservoir_simulator::get_main_loop ()
  {
    if (cm->n_phases == 3)
      return new main_loop_calc <true, true, true> (this);
    else if (cm->is_water () && cm->is_oil ())
      return new main_loop_calc <true, false, true> (this);
    else if (cm->is_gas () && cm->is_oil ())
      return new main_loop_calc <false, true, true> (this);
    else if (cm->is_water ())
      return new main_loop_calc <true, false, false> (this);
    else if (cm->is_gas ())
      return new main_loop_calc <false, true, false> (this);
    else if (cm->is_oil ())
      return new main_loop_calc <false, false, true> (this);
    else
      bs_throw_exception ("Unknown phase model");
  }

  void
  reservoir_simulator::main_loop ()
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

  void
  reservoir_simulator::pre_large_step (const sp_event_base_list_t &event_list)
  {
    mloop->apply_events (event_list);
    reservoir_->pre_large_step (cm, hdm_->mesh);
  }

  std::string
  reservoir_simulator::model_filename () const
  {
    return model_filename_;
  }

  // create object
  BLUE_SKY_TYPE_STD_CREATE (reservoir_simulator)
  BLUE_SKY_TYPE_STD_COPY (reservoir_simulator)

  // array map implementation
  BLUE_SKY_TYPE_IMPL (reservoir_simulator, bs_node, "reservoir_simulator", "Reservoir simulator", "Reservoir simulator");



  bool
  reservoir_simulator_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, reservoir_simulator::bs_type ());    BS_ASSERT (res);

    return res;
  }
}
