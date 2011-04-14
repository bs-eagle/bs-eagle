import bs
import numpy

#g = bs.lsolvers.gmres()
m = bs.mx.bcsr_matrix()
t = bs.mx.bcsr_matrix_tools ()

#n=150
#t.random_init(m,n,1,0.1,5)

fname="spe10.csr"
t.ascii_read_from_csr_format (m, fname)
n = m.get_n_rows ()

rhs = numpy.ones (n)
sol = numpy.zeros (n)

amg = bs.swift.amg_solver ()

p = amg.get_prop ()
p.set_i ("n_pre_smooth_iters_idx", 1)

#pmis2 = bs.swift.pmis2_coarse ()
#amg.set_coarser (0, pmis2)

amg.setup(m)
#amg.solve(m,rhs,sol)
#print sol