import sys
import bs
import numpy

hdm = bs.bs_bos_core_data_storage.hdm()

hdm.get_prop().add_property_i (0, "mesh", "mesh type")
hdm.get_prop().set_i ("mesh", 1)

#hdm.init("lalala.h5")
hdm.read_keyword_file(sys.argv[1])

#hdm.read_keyword_file('G:\\projects\\UfaSolver\\tests\\Frac\\Test2.ppp.files\\bos\\model.data')


#pool = hdm.get_pool()
#pool.py_set_pool_dims ([10, 10, 10])
#act = pool.get_i_data("ACTNUM")

mesh = hdm.get_mesh()

mesh.init_props(hdm)

mesh.init_ext_to_int()

mesh.check_data()

#jac = bs.mx.bcsr_matrix()
#flux_conn = bs.bs_mesh.flux_conn()
#bound = numpy.arange(1)
#mesh.build_jacobian_and_flux_connections(jac, flux_conn, bound)

#t = bs.mx.bcsr_matrix_tools ()
#t.ascii_write_to_csr_format (jac, "reg", 0)


#arrays = mesh.calc_element_tops()
#tops = arrays[0]
#tops = tops.reshape((-1,3))
#indexes = arrays[1]
#poro = arrays[2]

#arrays = mesh.calc_element_center()
#centers = arrays[0]
#centers = centers.reshape((-1,3))
#poro = arrays[1]
