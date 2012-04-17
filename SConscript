Import("*")

if build_kind == 'init':
    # create includes if not exists
    try: includes
    except NameError :
        includes = dict()
    else :
        pass

    includes["kernel"]                      = ["#kernel/include", "#kernel/include/python"]
    includes["bs_bos_core_base"]            = Dir("bs_bos_core_base/include")
    includes["bs_bos_core_data_storage"]    = Dir("bs_bos_core_data_storage/include")
    includes["bs_bos_core"]                 = Dir("bs_bos_core/include")
    includes["bs_cpr_prec"]                 = Dir("../bs_cpr_prec/include")
    includes["bs_csr_ilu_prec"]             = Dir("bs_csr_ilu_prec/include")
    includes["bs_lsolvers"]                 = Dir("bs_lsolvers/include")
    includes["bs_mesh"]                     = Dir("bs_mesh/include")
    includes["bs_mesh_src"]                 = Dir("bs_mesh/src")
    includes["bs_mtx"]                      = Dir("bs_mtx/include")
    includes["bs_pvt"]                      = Dir("bs_pvt/include")
    includes["bs_scal"]                     = Dir("bs_scal/include")
    includes["bs_swift"]                    = Dir("bs_swift/include")
    includes["common_types"]                = Dir("common_types/include")
    includes["common_alg"]                  = Dir("common_alg/include")
    includes["fluids"]                      = Dir("fluids/include")
    includes["sql_well"]                    = Dir("sql_well/include")
    includes["sql_well_src"]                = Dir("sql_well/src")
    includes["hdm_fluid"]                   = Dir("hdm_fluid/include")
    includes["bs_hdf5_storage"]             = Dir("hdf5_storage/src")
    includes["bs_pvt_src"]                  = Dir("bs_pvt/src")
    includes["bs_scal_src"]                 = Dir("bs_scal/src")
    #includes["common_types_src"]            = Dir("common_types/src")
    Export ("includes")

else :
    SConscript ("bs_bos_core_base/SConscript.bs")
    SConscript ("common_types/SConscript.bs")
    SConscript ("common_alg/SConscript.bs")
    SConscript ("fluids/SConscript.bs")
    SConscript ("sql_well/SConscript.bs")
    SConscript ("bs_mtx/SConscript.bs")
    SConscript ("bs_lsolvers/SConscript.bs")
    SConscript ("bs_swift/SConscript.bs")
#    SConscript ("hdm_fluid/SConscript.bs")
    SConscript ("bs_bos_core_data_storage/SConscript.bs")
    SConscript ("bs_mesh/SConscript.bs")
    SConscript ("bs_scal/SConscript.bs")
    SConscript ("bs_pvt/SConscript.bs")
    SConscript ("bs_bos_core/SConscript.bs")
    SConscript ("hdf5_storage/SConscript.bs")
