import bs
import numpy

hdm = bs.bs_bos_core_data_storage.hdm()

hdm.get_prop().add_property_i (0, "mesh", "mesh type")
hdm.get_prop().set_i ("mesh", 1)

hdm.init()
#hdm.read_keyword_file('D:\\projects\\UfaSolver\\tests\\bos\\model.data')
#hdm.read_keyword_file('G:\\projects\\UfaSolver\\tests\\Frac\\Test2.ppp.files\\bos\\model.data')
exit()
pool = hdm.get_pool()

mesh = hdm.get_mesh()

mesh.init_props(hdm)

#mesh.init_ext_to_int()

#jac = bs.mx.bcsr_matrix_dif()
#flux_conn = bs.bs_mesh.flux_conn_dif()
#bound = numpy.arange(1)
#mesh.build_jacobian_and_flux_connections(jac, flux_conn, bound)


#arrays = mesh.calc_element_tops()
#tops = arrays[0]
#tops = tops.reshape((-1,3))
#indexes = arrays[1]
#poro = arrays[2]

arrays = mesh.calc_element_center()
centers = arrays[0]
centers = centers.reshape((-1,3))
poro = arrays[1]
