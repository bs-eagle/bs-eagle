#ifndef BS_HDM_MODEL_IFACE_H
#define BS_HDM_MODEL_IFACE_H
/*!
\file hdm.h
\brief class hdm
\author Morozov Andrey
*/

namespace blue_sky {

  class keyword_manager_iface;
  class bos_reader_iface;
  class idata;
  class rs_mesh_iface;
  class h5_pool_iface;
  class prop_iface;
  class scal_3p_iface;
  class pvt_3p_iface;
  class init_model_iface;
  class event_manager_iface;
  class well_pool_iface;
  class equil_model_iface;

  class BS_API_PLUGIN hdm_iface: public objbase
    {
    public:

      typedef smart_ptr <idata, true>                       sp_idata ;
      typedef smart_ptr <bos_reader_iface, true>            sp_reader_t;
      typedef smart_ptr <keyword_manager_iface, true>       sp_km_t;
      typedef smart_ptr <rs_mesh_iface, true>               sp_mesh_iface_t;
      typedef smart_ptr <h5_pool_iface, true>               sp_h5_pool_t;
      typedef smart_ptr <prop_iface, true>                  sp_prop_t;

      //METHODS
      virtual ~hdm_iface() {};

    public:
    
      // initialize fluids
      virtual void init_fluids(t_int n_pvt_regions, t_int n_scal_regions) = 0;
      // initialize fluids
      virtual void init_equil(t_int n_equil_regions) = 0;
    
      // initialize data manager
      virtual void init(const std::string &model_name, bool use_memory = false) = 0;
      
      // read keyword file 
      virtual void read_keyword_file(const std::string filename) = 0;
      
      // ACCESS
      
      // input data accessor
      virtual sp_idata get_data () = 0;
      
      // reader accessor
      virtual sp_reader_t get_reader () = 0;
      
      // keyword manager accessor
      virtual sp_km_t get_keyword_manager () = 0;

      // mesh accessor
      virtual sp_mesh_iface_t get_mesh () = 0;

      virtual sp_h5_pool_t get_pool () = 0;

      virtual sp_prop_t get_prop () = 0;

      virtual t_double get_darcy_constant () = 0;

      virtual BS_SP (scal_3p_iface) get_scal () = 0;
      
      virtual BS_SP (pvt_3p_iface) get_pvt () = 0;

      virtual BS_SP (init_model_iface) get_init_model () = 0;

      virtual BS_SP (event_manager_iface) get_event_manager () = 0;
      
      virtual BS_SP (well_pool_iface) get_well_pool () = 0;
  
      virtual BS_SP (equil_model_iface) get_equil_model ()  = 0;
      // SET

      // mesh setter
      virtual void set_mesh (sp_mesh_iface_t mesh) = 0;

      virtual void set_init_model (BS_SP (init_model_iface) model) = 0;
  
      virtual void set_equil_model (BS_SP (equil_model_iface) model) = 0;
    
    };

} //ns bs

#endif //BS_HDM_MODEL_IFACE_H
