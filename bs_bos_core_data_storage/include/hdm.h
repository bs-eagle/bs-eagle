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
#include "read_class.h"
#include "scal_3p_iface.hpp"
#include "scal_dummy_iface.h"
#include "pvt_dummy_iface.h"
#include "init_model_iface.hpp"

namespace blue_sky {

  
  class BS_API_PLUGIN hdm: public hdm_iface
    {
    public:

      typedef smart_ptr <idata, true>										      sp_idata ;
      typedef smart_ptr <FRead, true>													sp_reader_t;
      typedef smart_ptr <keyword_manager_iface, true>					sp_km_t;
      typedef smart_ptr <rs_mesh_iface, true>		              sp_mesh_iface_t;

      //METHODS
      ~hdm();

    public:
    
      // initialize data manager
      void init(const std::string &model_name);
      
      // read keyword file 
      void read_keyword_file(const std::string filename);
      
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

      t_double get_darcy_constant () {return ph_const.darcy_constant;};

      BS_SP (scal_3p_iface) get_scal () { return scal_3p_; }
      
      BS_SP (scal_dummy_iface) get_scal_dummy () { return scal_dummy; }
      
      BS_SP (pvt_dummy_iface) get_pvt_dummy () { return pvt_dummy; }

      BS_SP (init_model_iface) get_init_model () { return init_model_; }

  
      // SET

      void set_mesh (sp_mesh_iface_t new_mesh) {mesh = new_mesh;};

      void set_init_model (BS_SP (init_model_iface) model) { init_model_ = model; }
      
      //CHECK
      void check_arrays_for_inactive_blocks () const;
      //void update_geometry() const;
      

    public:
      BLUE_SKY_TYPE_DECL (hdm)

    public:
      sp_idata                  data;												//!< data storage
      sp_reader_t               reader;                     //!< parser
      sp_km_t                   km;                         //!< keyword manager
      sp_mesh_iface_t           mesh;                       //!< mesh
      BS_SP (scal_3p_iface)     scal_3p_;
      BS_SP (scal_dummy_iface)  scal_dummy;
      BS_SP (pvt_dummy_iface)   pvt_dummy;
      BS_SP (init_model_iface)  init_model_;
      locale_keeper             lkeeper;
      physical_constants        ph_const;                //!< default physical constants
    };

} //ns bs
#endif //BS_HDM_MODEL_H
