/** 
 * @file py_sql_well.cpp
 * @brief python wraper for well fracture storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include "py_sql_well.h"
#include "sql_well.h"
#include "bs_prop_base.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  void py_export_compdat_ident();
  void py_export_sql_well_serialize();

  //////////////////////////////////////////////////////////////////////////

  std::string sqw_get_file_name(const sql_well& wp) {
    return wp.file_name;
  }
  void sqw_set_file_name(sql_well& wp, const std::string& new_fname) {
    wp.file_name = new_fname;
  }

  // old-fashioned overload for 3 params
  //int sqw_add_branch_gis_old(sql_well& wp, const std::string &wname, const std::string &branch,
  //                                sql_well::sp_gis_t g) {
  //  return wp.add_branch_gis(wname, branch, g);
  //}

  //sql_well::sp_gis_t sqw_get_branch_gis_old(sql_well& wp, const std::string &wname,
  //  const std::string &branch) {
  //  return wp.get_branch_gis(wname, branch);
  //}

  template< typename T >
  struct sql_well_exporter_plus {
    template< typename class_t >
    static class_t &
    export_class(class_t& class__) {
      py_sql_well_exporter< T >::export_class(class__)
        .add_property("file_name", &sqw_get_file_name, &sqw_set_file_name)
        //.def ("add_branch_gis", &sqw_add_branch_gis_old,
        //    args ("well_name", "branch_name", "gis"), "Add gis to the well branch")
        //.def ("get_branch_gis", &sqw_get_branch_gis_old,
        //    args ("well_name", "branch_name"), "Get gis for the well branch")
        ;
      return class__;
    }
  };

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_sql_well ()
  {
    using namespace boost::python;

    base_exporter <well_pool_iface, py_sql_well_exporter>::export_class ("well_pool_iface");

    class_exporter <sql_well, well_pool_iface, sql_well_exporter_plus>::export_class ("sql_well");

    py_export_compdat_ident();

    py_export_sql_well_serialize();

    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< sql_well >,
      smart_ptr< well_pool_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
