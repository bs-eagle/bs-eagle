/** 
 * @file py_sql_well.cpp
 * @brief python wraper for well fracture storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include "py_sql_well.h"
#include "sql_well.h"
#include "bs_serialize.h"
#include "bs_prop_base.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

BLUE_SKY_TYPE_SERIALIZE_GUID(blue_sky::sql_well)
BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::sql_well)

namespace blue_sky {
namespace python {
  void py_export_compdat_ident();

  //////////////////////////////////////////////////////////////////////////

  // helper function to write well_pool to given fname
  std::string serialize_to_str_fname(
      smart_ptr< well_pool_iface > wp,
      const std::string& prj_path,
      const std::string& prj_name
  ) {
    // save project path for serialization code
    kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(BS_KERNEL.find_type("hdm").td_);
    //std::string well_pool_filename = prj_name + "_well_pool.db";
    p_dt->insert< std::string >(prj_path);
    // and project name
    p_dt->insert< std::string >(prj_name);

    // invoke serializetion
    std::string res = serialize_to_str< well_pool_iface >(wp);
    // clear table
    p_dt->clear< std::string >();

    return res;
  }

  std::string sqw_get_file_name(const sql_well& wp) {
    return static_cast< const sql_well& >(wp).file_name;
  }
  void sqw_set_file_name(sql_well& wp, const std::string new_fname) {
    static_cast< sql_well& >(wp).file_name = new_fname;
  }

  template< typename T >
  struct sql_well_exporter_plus {
    template< typename class_t >
    static class_t &
    export_class(class_t& class__) {
      py_sql_well_exporter< T >::export_class(class__)
        .add_property("file_name", &sqw_get_file_name, &sqw_set_file_name);
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


    std::string (*s2s_wpi)(smart_ptr< well_pool_iface, true >&) = &blue_sky::serialize_to_str< well_pool_iface >;
    std::string (*s2s_sqw)(smart_ptr< sql_well, true >&) = &blue_sky::serialize_to_str< sql_well >;
    smart_ptr< well_pool_iface, true > (*sfs_wpi)(const std::string&) = &blue_sky::serialize_from_str< well_pool_iface >;
    smart_ptr< sql_well, true > (*sfs_sqw)(const std::string&) = &blue_sky::serialize_from_str< sql_well >;

    def("serialize_to_str", s2s_wpi);
    def("serialize_to_str", s2s_sqw);
    def("serialize_from_str", sfs_wpi);
    def("serialize_from_str", sfs_sqw);

    //def("serialize_to_str", &blue_sky::serialize_to_str< well_pool_iface >);
    //def("serialize_from_str", &blue_sky::serialize_from_str< well_pool_iface >);
    def("serialize_to_str_fname", &serialize_to_str_fname);

    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< sql_well >,
      smart_ptr< well_pool_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
