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
    includes["bs_mtx"]                      = Dir("bs_mtx/include")
    includes["bs_pvt"]                      = Dir("bs_pvt/include")
    includes["bs_scal"]                     = Dir("bs_scal/include")
    includes["bs_swift"]                    = Dir("bs_swift/include")
    includes["common_types"]                = Dir("common_types/include")
    includes["hdm_fluid"]                   = Dir("hdm_fluid/include")
    Export ("includes")

else :
    SConscript ("bs_bos_core_base/SConscript.bs")
    SConscript ("common_types/SConscript.bs")
    SConscript ("bs_mtx/SConscript.bs")
    SConscript ("bs_lsolvers/SConscript.bs")
    SConscript ("bs_swift/SConscript.bs")
    SConscript ("hdm_fluid/SConscript.bs")
    SConscript ("bs_bos_core_data_storage/SConscript.bs")
    SConscript ("bs_mesh/SConscript.bs")
    SConscript ("bs_scal/SConscript.bs")
    SConscript ("bs_pvt/SConscript.bs")
    SConscript ("bs_bos_core/SConscript.bs")
    SConscript ("std_well/SConscript.bs")
