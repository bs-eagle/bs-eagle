#ifndef BS_WELL_STORAGE_H
#define BS_WELL_STORAGE_H
/**
 * \file well_storage.h
 * \brief multi-index well storage 
 * \author Mark Khait
 * \date 07.07.2011
 * */
 
 
#include "well_storage_iface.h"
#include "well_iface.h"

namespace blue_sky {

  class BS_API_PLUGIN well_storage: public well_storage_iface
    {
    public:
    
      typedef BS_SP (well_iface) sp_well_iface; 
  
      //METHODS
      ~well_storage();

    public:
    
      //typedef std::map <std::string, std::list<std::string>> str_map;
      
      void add_group (group_params *params);
      void add_well (well_params *params, const std::string &group);
      void add_conng (conng_params *params, const std::string &well);
      
      int add_well_to_group (const std::string &well, const std::string &group);
      
      well_params* get_well_params (const std::string &well, t_double date);
      group_params* get_group_params (const std::string &group, t_double date);
      
      void set_group_fparam (const std::string &group, t_double date, int col, t_double val);
      void set_group_iparam (const std::string &group, t_double date, int col, t_int val);
      
      void set_group_fparams (const std::string &group, t_double date, spv_double);
      void set_group_iparams (const std::string &group, t_double date, spv_int);
      
      void set_well_fparam (const std::string &well, t_double date, int col, t_double val);
      void set_well_iparam (const std::string &well, t_double date, int col, t_int val);
      
      void set_well_fparams (const std::string &well, t_double date, spv_double);
      void set_well_iparams (const std::string &well, t_double date, spv_int);
      
      //void update_group (const std::string &group);
      
      template <typename data_set, int ni, int nf>
      spv_double get_dates (const std::string &name, data_set data);
      
      template <typename data_set, int ni, int nf>
       std::list<std::string> get_names (data_set data);
      
      std::list<std::string> get_group_names ();
      spv_double get_group_dates (const std::string &group);
      spv_double get_group_fvalues (const std::string &group, int col);
      spv_int get_group_ivalues (const std::string &group, int col);
      
      std::list<std::string> get_well_names ();
      spv_double get_well_dates (const std::string &well);
      spv_double get_well_fvalues (const std::string &well, int col);
      spv_int get_well_ivalues (const std::string &well, int col);
      
      
      
       /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const 
      { 
        return "To be continued...\n";
      };
      
    private:
      int add_child (const std::string &child, const std::string &parent, str_map &parent_map);
      
      template <typename params_t, typename data_set_t, int ni, int nf>
        params_t* get_params (data_set_t &data, const std::string &well, t_double date);
      
      template <typename data_set_t, typename vector_t, int ni, int nf>
        blue_sky::smart_ptr< vector_t > get_values (const std::string &name, data_set_t &data, int col, int type);
        
        
      //! set exact param on specified date
      //! data - target multiindex
      //! name - target object name
      //! date - target date
      //! col - number of target param
      //! val - new value (operation = 0) or value to add (operation = 1)
      //! type - parameter type, floating point with type = 0 and integer with type = 1
      //! operation - operation type, assignment with operation = 0 and addition with operation = 1
      template <typename params_t, typename data_t, typename data_set_t, int ni, int nf>
        void set_param (data_set_t &data, const std::string &name, t_double date, int col, data_t &val, int type, int operation);
      
      template <typename params_t, typename vector_t, typename data_set_t, int ni, int nf>
        void set_params (data_set_t &data, const std::string &name, t_double date, blue_sky::smart_ptr< vector_t > vals, int type, int operation);
        
      //! set exact param to child and take this change into account for all parents
      template <typename child_params_t, typename parent_params_t, typename data_t,typename child_data_set_t, typename parent_data_set_t, 
                int child_ni, int child_nf, int parent_ni, int parent_nf>
        void set_child_param (child_data_set_t &child_data, parent_data_set_t &parent_data, str_map &parent_map, const std::string &child_name, 
                                           t_double date, int col, data_t val, int type);
        
      template <typename child_params_t, typename parent_params_t, typename vector_t,typename child_data_set_t, typename parent_data_set_t, 
                int child_ni, int child_nf, int parent_ni, int parent_nf>
        void set_child_params (child_data_set_t &child_data, parent_data_set_t &parent_data, str_map &parent_map, const std::string &child_name, 
                             t_double date, blue_sky::smart_ptr< vector_t > vals, int type);
     
     
      
    public:
      BLUE_SKY_TYPE_DECL (well_storage)

    public:
    
      group_set groups;
      well_set wells;
      conng_set conngs;
      
      str_map groups_map;
      str_map wells_map;
      
      std::map<std::string, sp_well_iface> well_branches;
      
    };

} //ns bs
#endif //BS_WELL_STORAGE_H
