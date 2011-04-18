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
p.set_i ("n_last_level_points", 10)

pmis2 = bs.swift.pmis2_coarse ()
amg.set_coarser (0, pmis2)
pbuild2 = bs.swift.standart2_pbuild ()
amg.set_pbuilder (0, pbuild2)

#amg as solver
amg.setup(m)
#amg.solve(m,rhs,sol)
#print amg.get_prop()
m_lev = amg.get_matrices ()
cf_lev = amg.get_cf_markers ()
s_lev = amg.get_s_markers ()
p_lev = amg.get_p_matrices ()

for i in xrange(0, len (m_lev)):
    fname = "A_" + repr(i) + ".csr"
    t.ascii_write_to_csr_format (m_lev[i], fname, 0)
    
for i in xrange(0, len (p_lev)):
    fname = "P_" + repr(i) + ".csr"
    t.ascii_write_to_csr_format (p_lev[i], fname, 0)
    
for i in xrange(0, len (s_lev)):
    print "S_MARKERS LEVEL " + repr(i) + ":"
    print s_lev[i]

for i in xrange(0, len (cf_lev)):
    print "CF_MARKERS LEVEL " + repr(i) + ":"
    print cf_lev[i]
    

#amg as prec
#p.set_i ("maxiters", 1)
#g.set_prec(amg)
#g.setup(m)
#g.solve(m,rhs,sol)
#print g.get_prop()

#print sol