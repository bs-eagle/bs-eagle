Import("*")

if build_kind == 'init':
	# check that includes exists
	try : includes
	except NameError :
		includes = dict()
	else :
		pass

	# fill includes
	includes["kernel"]                      = ["#kernel/include", "#kernel/include/python"]
	includes["bs_bos_core_base"]            = Dir("bs_bos_core_base/include")
	includes["bs_bos_core_data_storage"]    = Dir("bs_bos_core_data_storage/include")
	includes["bs_bos_core"]                 = Dir("bs_bos_core/include")
	includes["bs_base_linear_solvers"]      = Dir("bs_base_linear_solvers/include")
	includes["bs_csr_ilu_prec"]             = Dir("bs_csr_ilu_prec/include")
	includes["bs_mesh"]                     = Dir("bs_mesh/include")
	includes["bs_matrix"]                   = Dir("bs_matrix/include")
	includes["bs_mtx"]                      = Dir("bs_mtx/include")
	includes["bs_pvt"]                      = Dir("bs_pvt/include")
	includes["bs_scal"]                     = Dir("bs_scal/include")
	Export ("includes")

else :
	SConscript ("bs_bos_core_base/SConscript.bs")
	SConscript ("bs_mtx/SConscript.bs")
	SConscript ("bs_matrix/SConscript.bs")
	SConscript ("bs_base_linear_solvers/SConscript.bs")
	SConscript ("bs_csr_ilu_prec/SConscript.bs")
	SConscript ("bs_bos_core_data_storage/SConscript.bs")
	SConscript ("bs_mesh/SConscript.bs")
	SConscript ("bs_scal/SConscript.bs")
	SConscript ("bs_pvt/SConscript.bs")
	SConscript ("bs_bos_core/SConscript.bs")

