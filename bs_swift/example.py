import bs
import numpy

g = bs.lsolvers.gmres()
m = bs.mx.bcsr_matrix()
t = bs.mx.bcsr_matrix_tools ()

#n=150
#t.random_init (m, n, 1, 0.1, 5)

n=5
t.gen_2d_laplas (m, n)
n = m.get_n_rows ()

#fname="spe10.csr"
#t.ascii_read_from_csr_format (m, fname)
#n = m.get_n_rows ()

rhs = numpy.ones (n)
sol = numpy.zeros (n)

amg = bs.swift.amg_solver ()

p = amg.get_prop ()
p.set_i ("maxiters", 20)
p.set_i ("n_last_level_points", 15)

pmis2 = bs.swift.pmis2_coarse ()
amg.set_coarser (0, pmis2)

#amg.setup(m)
#amg.solve(m,rhs,sol)
#print amg.get_prop()

p.set_i ("maxiters", 1)
g.set_prec(amg)
g.setup(m)
g.solve(m,rhs,sol)
print g.get_prop()

#print sol