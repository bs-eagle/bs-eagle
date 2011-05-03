#ifndef BS_HDM_MODEL_IFACE_H
#define BS_HDM_MODEL_IFACE_H
/*!
\file hdm.h
\brief class hdm
\author Morozov Andrey
*/

namespace blue_sky {

  class keyword_manager_iface;
  class FRead;
  class idata;
  class rs_mesh_iface;
  class h5_pool_iface;
  class prop_iface;

  class BS_API_PLUGIN hdm_iface: public objbase
    {
    public:

      typedef smart_ptr <idata, true>							      sp_idata ;
      typedef smart_ptr <FRead, true>								    sp_reader_t;
      typedef smart_ptr <keyword_manager_iface, true>		sp_km_t;
      typedef smart_ptr <rs_mesh_iface, true>		        sp_mesh_iface_t;
      typedef smart_ptr <h5_pool_iface, true>						sp_h5_pool_t;
      typedef smart_ptr <prop_iface, true>							sp_prop_t;

      //METHODS
      virtual ~hdm_iface() {};

    public:
    
      // initialize data manager
      virtual void init() = 0;
      
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

  
      // SET

      // mesh setter
      virtual void set_mesh (sp_mesh_iface_t mesh) = 0;
      
    };

} //ns bs

#endif //BS_HDM_MODEL_IFACE_H
