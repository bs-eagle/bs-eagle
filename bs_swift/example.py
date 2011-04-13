import bs
import numpy

#g = bs.lsolvers.gmres()
m = bs.mx.bcsr_matrix()
t = bs.mx.bcsr_matrix_tools ()

#fname=spe10.csr
#t.ascii_read_from_csr_format (m, fname)
n=50
t.random_init(m,n,1,0.1,5)
rhs = numpy.ones (n)
sol = numpy.zeros (n)

amg = bs.swift.amg_solver()
amg.setup(m)
#amg.solve(m,rhs,sol)
#print sol