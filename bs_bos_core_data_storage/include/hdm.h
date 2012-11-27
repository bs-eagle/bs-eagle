#ifndef BS_HDM_MODEL_H
#define BS_HDM_MODEL_H
/*!
\file hdm.h
\brief class hdm
\author Morozov Andrey
*/

#include "data_class.h"
#include "keyword_manager.h"
#include "hdm_iface.h"
#include "locale_keeper.h"
#include "rs_mesh_iface.h"
#include "bos_reader_iface.h"
#include "scal_3p_iface.hpp"
#include "pvt_3p_iface.h"
#include "init_model_iface.hpp"
#include "event_manager_iface.hpp"
#include "well_pool_iface.h"
#include "equil_model_iface.h"

namespace blue_sky {

  
  class BS_API_PLUGIN hdm: public hdm_iface
    {
    public:

      typedef smart_ptr <idata, true>                               sp_idata ;
      typedef smart_ptr <bos_reader_iface, true>                    sp_reader_t;
      typedef smart_ptr <keyword_manager_iface, true>               sp_km_t;
      typedef smart_ptr <rs_mesh_iface, true>                       sp_mesh_iface_t;

      //METHODS
      ~hdm();

    public:
    
      // initialize fluids
      void init_fluids(t_int n_pvt_regions, t_int n_scal_regions);
      void init_equil (t_int n_equil_regions);
      
      // initialize data manager
      void init(const std::string &model_name);

      // initialize process params
      void init_proc_params();
      
      // read keyword file 
      void read_keyword_file(const std::string filename);
      
      void test_well_storage ();
      
      // ACCESS
      
      // keyword manager accessor
      sp_idata get_data () {return data;};
      
      // keyword manager accessor
      sp_reader_t get_reader () {return reader;};
      
      // keyword manager accessor
      sp_km_t get_keyword_manager () {return km;};

      // mesh accessor
      sp_mesh_iface_t get_mesh () {return mesh;};

      sp_h5_pool_t get_pool () {return data->h5_pool;};

      sp_prop_t get_prop () {return data->props;};

      sp_prop_t get_proc_params () {return data->proc_params;};

      t_double get_darcy_constant () {return ph_const.darcy_constant;};

      BS_SP (scal_3p_iface) get_scal () { return scal_3p_; }
      
      BS_SP (pvt_3p_iface) get_pvt () { return pvt_3p_; }

      BS_SP (init_model_iface) get_init_model () { return init_model_; }

      BS_SP (event_manager_iface) get_event_manager () { return event_manager_; }
      
      BS_SP (well_pool_iface) get_well_pool () { return well_pool_; }
      
      BS_SP (equil_model_iface) get_equil_model ()  { return equil_model_; }
      // SET

      void set_mesh (sp_mesh_iface_t new_mesh) {mesh = new_mesh;};

      void set_init_model (BS_SP (init_model_iface) model) { init_model_ = model; }

      void set_equil_model (BS_SP (equil_model_iface) model) { equil_model_ = model; }
      
      //CHECK
      void check_arrays_for_inactive_blocks () const;
      //void update_geometry() const;
      

    public:
      BLUE_SKY_TYPE_DECL (hdm)

    public:
      sp_idata                  data;                                               //!< data storage
      sp_reader_t               reader;                     //!< parser
      sp_km_t                   km;                         //!< keyword manager
      sp_mesh_iface_t           mesh;                       //!< mesh
      BS_SP (scal_3p_iface)     scal_3p_;
      BS_SP (pvt_3p_iface)      pvt_3p_;
      BS_SP (equil_model_iface) equil_model_;
      BS_SP (init_model_iface)  init_model_;
      BS_SP (event_manager_iface) event_manager_;
      BS_SP (well_pool_iface)   well_pool_;                 //!< SQL well pool
      locale_keeper             lkeeper;
      physical_constants        ph_const;                //!< default physical constants
    };

} //ns bs
#endif //BS_HDM_MODEL_H
