#ifndef BS_WELL_STORAGE_IFACE_H
#define BS_WELL_STORAGE_IFACE_H
/**
 * \file well_storage_iface.h
 * \brief interface for multi-index well storage 
 * \author Mark Khait
 * \date 12.07.2011
 * */

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
    
#include "conf.h"

using namespace boost;
using namespace boost::multi_index;    
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace std; 



namespace blue_sky {


  #define NI_GROUP 10
  #define NF_GROUP 10
  #define NI_WELL 20
  #define NF_WELL 20
  #define NI_CONNG 15
  #define NF_CONNG 15
  
  template<int n_i_params, int n_f_params> class params
    {
      public:
      
        std::string name;
        t_double      date;
        
        int         i_params[n_i_params];
        t_double    f_params[n_f_params];
    };
    
  typedef params<NI_GROUP,NF_GROUP> group_params;
  typedef params<NI_WELL,NF_WELL> well_params;
  typedef params<NI_CONNG,NF_CONNG> conng_params;
  
  template<int n_i_params, int n_f_params> class comp_date
  {
    public:  
      bool operator()(const params<n_i_params,n_f_params>& w1,const params<n_i_params,n_f_params>& w2)const 
      {
        if (w1.name < w2.name)
         return true;
        else if (w1.name == w2.name)
          return (w1.date < w2.date);
        else
          return false;
      } // && w1.date < w2.date
    };
     
  template<int n_i_params, int n_f_params> class comp_name
      {
        public:  
          bool operator()(std::string x,const params<n_i_params,n_f_params>& w2)const{return x < w2.name;}
          bool operator()(const params<n_i_params,n_f_params>& w1,std::string x)const{return w1.name < x;}
      };

  typedef multi_index_container<
    group_params,
    indexed_by<
      ordered_unique<identity<group_params>, comp_date<NI_GROUP,NF_GROUP> >,
      ordered_non_unique<member<group_params, t_double, &group_params::date> >
    > 
  > group_set;
  
  typedef multi_index_container<
    well_params,
    indexed_by<
      ordered_unique<identity<well_params>, comp_date<NI_WELL,NF_WELL> >,
      ordered_non_unique<member<well_params, t_double, &well_params::date> >
    > 
  > well_set;
  
  typedef multi_index_container<
    conng_params,
    indexed_by<
      ordered_unique<identity<conng_params>, comp_date<NI_CONNG,NF_CONNG> >,
      ordered_non_unique<member<conng_params, t_double, &conng_params::date> >
    > 
  > conng_set;
  
    
  class BS_API_PLUGIN well_storage_iface: public objbase
    {
    public:

      //METHODS
      virtual ~well_storage_iface() {};

    public:
    
      typedef std::map <std::string, std::set<std::string> > str_map;
      
      virtual void add_group (group_params *params) = 0;
      virtual void add_well (well_params *params, const std::string &group) = 0;
      virtual void add_conng (conng_params *params, const std::string &well) = 0;
      
      virtual int add_well_to_group (const std::string &well, const std::string &group) = 0;
      
      virtual void set_group_fparam (const std::string &group, t_double date, int col, t_double val) = 0;
      virtual void set_group_iparam (const std::string &group, t_double date, int col, t_int val) = 0;
      
      virtual void set_group_fparams (const std::string &group, t_double date, spv_double vals) = 0;
      virtual void set_group_iparams (const std::string &group, t_double date, spv_int vals) = 0;
      
      virtual void set_well_fparam (const std::string &well, t_double date, int col, t_double val) = 0;
      virtual void set_well_iparam (const std::string &well, t_double date, int col, t_int val) = 0;
      
      virtual void set_well_fparams (const std::string &well, t_double date, spv_double vals) = 0;
      virtual void set_well_iparams (const std::string &well, t_double date, spv_int vals) = 0;
      
      //virtual void update_group (const std::string &group) = 0;
      
      virtual spv_double get_group_dates (const std::string &group) = 0;
      virtual spv_double get_group_fvalues (const std::string &group, int col) = 0;
      virtual spv_int get_group_ivalues (const std::string &group, int col) = 0;
      
      virtual spv_double get_well_dates (const std::string &well) = 0;
      virtual spv_double get_well_fvalues (const std::string &well, int col) = 0;
      virtual spv_int get_well_ivalues (const std::string &well, int col) = 0;
      
      
    };

} //ns bs

#endif //BS_WELL_STORAGE_H
