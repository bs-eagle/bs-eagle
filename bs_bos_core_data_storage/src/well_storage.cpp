#include "bs_bos_core_data_storage_stdafx.h"

#include "well_storage.h"
//#include "main_def.h"

#include "write_time_to_log.h"


 
namespace blue_sky
  {
  
  well_storage::~well_storage ()
  {

  }

  
  well_storage::well_storage(bs_type_ctor_param /*param*/)
  {
  }

  well_storage::well_storage(const well_storage& src)
    :bs_refcounter ()
  {
    *this = src;
  }
  
  void well_storage::add_group (group_params *params)
  {
    groups.insert (*params);
  }
  
  void well_storage::add_well (well_params *params, const std::string &group)
  {
    str_map::iterator it;
    wells.insert (*params);
    
    if (group != "")
      {
        add_child (params->name, group, groups_map);
      }  
  }
  
  void well_storage::add_conng (conng_params *params, const std::string &well)
  {
    str_map::iterator it;
    std::set<std::string> str_set;
    
    conngs.insert (*params);
    if (well != "")
      {
         add_child (params->name, well, wells_map);
      }
  }
  
  int well_storage::add_well_to_group (const std::string &well, const std::string &group)
  {
    return add_child (well, group, groups_map);
  }
  
  int well_storage::add_child (const std::string &child, const std::string &parent, str_map &parent_map)
  {
    str_map::iterator it;
    std::set<std::string> str_set;
    
    it = parent_map.find(parent);
    if (it == parent_map.end ())
      {
        str_set.insert(child);
        return parent_map.insert (str_map::value_type(parent, str_set)).second;
      }
    else
      {
        it->second.insert(child);
        return 1;
      }
  }
  
  template <typename data_set, int ni, int nf>
  std::list<std::string> well_storage::get_names (data_set data)
  {
    typename data_set::template nth_index<0>::type& name_index = data.get<0>();
    typename data_set::template nth_index<0>::type::iterator name_it, name_b, name_e;
    std::list<std::string> names;
    long n = 0;
     
    name_b = name_index.begin();
    name_e = name_index.end();
    n = distance(name_b, name_e);
    
    if (n)
      {
        names.resize (n);
       
        for (name_it = name_b; name_it != name_e; ++name_it)
          {
            names.push_back(name_it->name);
          }
      }
    return names;
  }
   
  template <typename data_set, int ni, int nf>
  spv_double well_storage::get_dates (const std::string &name, data_set data)
  {
    typename data_set::template nth_index<0>::type& name_index = data.get<0>();
    typename data_set::template nth_index<0>::type::iterator name_it, name_b, name_e;
    spv_double res;
    t_double *res_data;
    long n = 0;
     
    name_b = name_index.lower_bound (name, comp_name<ni,nf>());
    name_e = name_index.upper_bound (name, comp_name<ni,nf>());
    n = distance(name_b, name_e);
    
    res = give_kernel::Instance().create_object(v_double::bs_type());
    if (n)
      {
        res->resize (n);
        res_data = &(*res)[0];
        n = 0;
        
        for (name_it = name_b; name_it != name_e; ++name_it)
          {
            res_data[n] = name_it->date;
            n++;
          }
      }
    return res;
  }
  
  template <typename data_set_t, typename vector_t, int ni, int nf>
  blue_sky::smart_ptr< vector_t > well_storage::get_values (const std::string &name, data_set_t &data, int col, int type = 0)
  {
    typename data_set_t::template nth_index<0>::type& name_index = data.get<0>();
    typename data_set_t::template nth_index<0>::type::iterator name_it, name_b, name_e;
    blue_sky::smart_ptr< vector_t > res;
    typename vector_t::value_type *res_data;
    long n = 0;
     
    name_b = name_index.lower_bound (name, comp_name<ni,nf>());
    name_e = name_index.upper_bound (name, comp_name<ni,nf>());
    n = distance(name_b, name_e);
    
    res = give_kernel::Instance().create_object(vector_t::bs_type());
    if (n)
      {
        res->resize (n);
        res_data = &(*res)[0];
        n = 0;
        
        if (type == 0)
          for (name_it = name_b; name_it != name_e; ++name_it)
            {
              res_data[n] = name_it->f_params[col];
              n++;
            }
        else
          for (name_it = name_b; name_it != name_e; ++name_it)
            {
              res_data[n] = name_it->i_params[col];
              n++;
            }
      }  
    return res;
  }
  
  
  template <typename params_t, typename data_set_t, int ni, int nf>
  params_t*
  well_storage::get_params (data_set_t &data, const std::string &well, t_double date)
  {
    typename data_set_t::template nth_index<0>::type& name_index = data.get<0>();
    typename data_set_t::iterator name_it, name_b, name_e;
     
    name_b = name_index.lower_bound (well, comp_name<ni,nf>());
    name_e = name_index.upper_bound (well, comp_name<ni,nf>());
    
    for (name_it = name_b; name_it != name_e; ++name_it)
      {
        if (name_it->date == date)
          {
            return &const_cast<params_t&>(*name_it);
          }
      }
    return 0;
  }
  
  
  template <typename params_t, typename data_t, typename data_set_t, int ni, int nf>
  void well_storage::set_param (data_set_t &data, const std::string &name, t_double date, int col, data_t &val, int type = 0, int operation = 0)
  {
    params_t *params;
    data_t tmp;
     
    params = get_params<params_t, data_set_t, ni, nf>(data, name, date);
    
    if (params)
      {
        if (type == 0)
          if (operation == 0)
            {
              tmp = params->f_params[col];
              params->f_params[col] = val;
              val -= tmp;
            }
          else
            params->f_params[col] += val;
        else
          if (operation == 0)
            {
              tmp = params->i_params[col];
              params->i_params[col] = val;
              val -= tmp;
            }
          else
            params->i_params[col] += val;
      }
    else
      {
        params = new params_t;
        params->name = name;
        params->date = date;
        if (type == 0)
          params->f_params[col] = val;
        else
          params->i_params[col] = val;
        data.insert (*params);
      }
  }
  
  template <typename params_t, typename vector_t, typename data_set_t, int ni, int nf>
  void well_storage::set_params (data_set_t &data, const std::string &name, t_double date, blue_sky::smart_ptr< vector_t > vals, int type = 0, int operation = 0)
  {
    params_t *params;
    typename vector_t::value_type *vals_data;
    typename vector_t::value_type tmp;
    int i;
     
    params = get_params<params_t, data_set_t, ni, nf>(data, name, date);
    vals_data = &(*vals)[0];
    
    if (params)
      {
        
        if (type == 0)
          if (operation == 0)
            for (i = 0; i < nf; i++)
              {
                tmp = params->f_params[i];
                params->f_params[i] = vals_data[i];
                vals_data[i] -= tmp;
              }
          else
            for (i = 0; i < nf; i++)
              params->f_params[i] += vals_data[i];
        else
          if (operation == 0)
            for (i = 0; i < ni; i++)
              {
                tmp = params->i_params[i];
                params->i_params[i] = vals_data[i];
                vals_data[i] -= tmp;
              }
          else
            for (i = 0; i < ni; i++)
            params->i_params[i] += vals_data[i];
      }
    else
      {
        params = new params_t;
        params->name = name;
        params->date = date;
        if (type == 0)
          for (i = 0; i < nf; i++)
            params->f_params[i] = vals_data[i];
        else
          for (i = 0; i < ni; i++)
            params->i_params[i] = vals_data[i];
        data.insert (*params);
      }
  }
  
  template <typename child_params_t, typename parent_params_t, typename data_t,typename child_data_set_t, typename parent_data_set_t, 
            int child_ni, int child_nf, int parent_ni, int parent_nf>
  void well_storage::set_child_param (child_data_set_t &child_data, parent_data_set_t &parent_data, str_map &parent_map, const std::string &child_name, 
                                           t_double date, int col, data_t val, int type = 0)
  {
    str_map::iterator it, b, e;
    
    // set child params and get the difference (inside vals)
    set_param <child_params_t, data_t, child_data_set_t, child_ni, child_nf> (child_data, child_name, date, col, val, type);
    
    // add the difference to all parents
    
    b = parent_map.begin();
    e = parent_map.end();
    for (it = b; it != e; ++it)
      {
        if (it->second.find(child_name) != it->second.end())
          set_param <parent_params_t, data_t, parent_data_set_t, parent_ni, parent_nf> (parent_data, it->first, date, col, val, type, 1); 
      }
  }
  
  template <typename child_params_t, typename parent_params_t, typename vector_t,typename child_data_set_t, typename parent_data_set_t, 
            int child_ni, int child_nf, int parent_ni, int parent_nf>
  void well_storage::set_child_params (child_data_set_t &child_data, parent_data_set_t &parent_data, str_map &parent_map, const std::string &child_name, 
                                           t_double date, blue_sky::smart_ptr< vector_t > vals, int type = 0)
  {
    str_map::iterator it, b, e;
    
    // set child params and get the difference (inside vals)
    set_params <child_params_t, vector_t, child_data_set_t, child_ni, child_nf> (child_data, child_name, date, vals, type);
    
    // add the difference to all parents
    
    b = parent_map.begin();
    e = parent_map.end();
    for (it = b; it != e; ++it)
      {
        if (it->second.find(child_name) != it->second.end())
          set_params <parent_params_t, vector_t, parent_data_set_t, parent_ni, parent_nf> (parent_data, it->first, date, vals, type, 1); 
      }
  }
  
  
  
  
  void well_storage::set_well_fparam (const std::string &well_name, t_double date, int col, t_double val)
  {
    set_child_param <well_params, group_params, t_double, well_set, group_set, NI_WELL, NF_WELL, NI_GROUP, NF_GROUP> (wells, groups, groups_map, well_name, date, col, val);
    well_branches.insert(std::pair<std::string, sp_well_iface> (well_name, BS_KERNEL.create_object ("well")));
  }
  
    
  void well_storage::set_well_iparam (const std::string &well_name, t_double date, int col, t_int val)
  {
    set_child_param <well_params, group_params, t_int, well_set, group_set, NI_WELL, NF_WELL, NI_GROUP, NF_GROUP> (wells, groups, groups_map, well_name, date, col, val, 1);
    well_branches.insert(std::pair<std::string, sp_well_iface> (well_name, BS_KERNEL.create_object ("well")));
  }
   
  void well_storage::set_well_fparams (const std::string &well_name, t_double date, spv_double vals)
  {
    set_child_params <well_params, group_params, v_double, well_set, group_set, NI_WELL, NF_WELL, NI_GROUP, NF_GROUP> (wells, groups, groups_map, well_name, date, vals);
    well_branches.insert(std::pair<std::string, sp_well_iface> (well_name, BS_KERNEL.create_object ("well")));
  }
  
  
  void well_storage::set_well_iparams (const std::string &well_name, t_double date, spv_int vals)
  {
    set_child_params <well_params, group_params, v_int, well_set, group_set, NI_WELL, NF_WELL, NI_GROUP, NF_GROUP> (wells, groups, groups_map, well_name, date, vals, 1);
    well_branches.insert(std::pair<std::string, sp_well_iface> (well_name, BS_KERNEL.create_object ("well")));
  }
  
  /*
  void well_storage::update_group (const std::string &group)
  {
    str_map::iterator gr;
    std::set<std::string>::iterator it, b, e;
    spv_double dates;
    
    gr = groups_map.find (group)
    b = gr->second.begin();
    e = gr->second.end();
    for (it = b; it != e; ++it)
      {
        dates = get_well_dates
      }
  }
  */
  
  void well_storage::set_group_fparam (const std::string &group_name, t_double date, int col, t_double val)
  {
    set_param <group_params, t_double, group_set, NI_GROUP, NF_GROUP> (groups, group_name, date, col, val);
  }
  
    
  void well_storage::set_group_iparam (const std::string &group_name, t_double date, int col, t_int val)
  {
    set_param <group_params, t_int, group_set, NI_GROUP, NF_GROUP> (groups, group_name, date, col, val, 1);
  }
   
  
  void well_storage::set_group_fparams (const std::string &group_name, t_double date, spv_double vals)
  {
    set_params <group_params, v_double, group_set, NI_GROUP, NF_GROUP> (groups, group_name, date, vals);
  }
  
  
  void well_storage::set_group_iparams (const std::string &group_name, t_double date, spv_int vals)
  {
    set_params <group_params, v_int, group_set, NI_GROUP, NF_GROUP> (groups, group_name, date, vals, 1);
  }
  
  std::list<std::string> well_storage::get_group_names ()
  {
    return get_names<group_set, NI_GROUP, NF_GROUP> (groups);
  }
  
    spv_double well_storage::get_group_dates (const std::string &group_name)
  {
    return get_dates<group_set, NI_GROUP, NF_GROUP> (group_name, groups);
  }
  
  spv_double well_storage::get_group_fvalues (const std::string &group_name, int col)
  {
    return get_values<group_set, v_double, NI_GROUP, NF_GROUP> (group_name, groups, col);
  }
  
  spv_int well_storage::get_group_ivalues (const std::string &group_name, int col)
  {
    return get_values<group_set, v_int, NI_GROUP, NF_GROUP> (group_name, groups, col, 1);
  }
  
  std::list<std::string> well_storage::get_well_names ()
  {
    return get_names<well_set, NI_WELL, NF_WELL> (wells);
  }
  
  spv_double well_storage::get_well_dates (const std::string &well_name)
  {
    return get_dates<well_set, NI_WELL, NF_WELL> (well_name, wells);
  }
  
  spv_double well_storage::get_well_fvalues (const std::string &well_name, int col)
  {
    return get_values<well_set, v_double, NI_WELL, NF_WELL> (well_name, wells, col);
  }
  
  spv_int well_storage::get_well_ivalues (const std::string &well_name, int col)
  {
    return get_values<well_set, v_int, NI_WELL, NF_WELL> (well_name, wells, col, 1);
  }

  
  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(well_storage)
  BLUE_SKY_TYPE_STD_COPY(well_storage)

  BLUE_SKY_TYPE_IMPL(well_storage, objbase, "well_storage", "BOS_Core well storage", "BOS_Core well storage")
}//ns bs
