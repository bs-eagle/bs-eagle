import bs
import numpy

g = bs.lsolvers.gmres()
m = bs.mx.bcsr_matrix()
t = bs.mx.bcsr_matrix_tools ()

#n=150
#t.random_init (m, n, 1, 0.1, 5)

n=50
t.gen_2d_laplas (m, n)
n = m.get_n_rows ()
rhs = numpy.ones (n)
sol = numpy.zeros (n)

# read system from file
#fname="spe10.csr"
#t.ascii_read_from_csr_format (m, fname)
#n = m.get_n_rows ()
#sol = numpy.zeros (n)
#fname="spe10.rhs"
#rhs = numpy.loadtxt (fname)

amg = bs.swift.amg_solver ()

p = amg.get_prop ()
p.set_f ("tolerance", 1e-5)
p.set_i ("maxiters", 25)
p.set_i ("n_last_level_points", 100)

pmis2 = bs.swift.pmis2_coarse ()
amg.set_coarser (0, pmis2)
pbuild2 = bs.swift.standart2_pbuild ()
amg.set_pbuilder (0, pbuild2)

# amg as solver
amg.setup (m)
amg.solve (m, rhs, sol)
print amg.get_prop()

n=60
t.gen_2d_laplas (m, n)
n = m.get_n_rows ()
rhs = numpy.ones (n)
sol = numpy.zeros (n)
amg.setup (m)
amg.solve (m, rhs, sol)
print amg.get_prop()


exit ()

# amg as prec
gp = g.get_prop ()
gp.set_i ("maxiters", 15)
p.set_i ("maxiters", 1)
g.set_prec (amg)
g.setup (m)
g.solve (m, rhs, sol)
print g.get_prop ()
#print amg.get_prop ()

exit()

# output AMG solve data:
sol_lev = amg.get_sol ()
for i in xrange(0, len (sol_lev)):
    fname = "sol." + repr(i)
    numpy.savetxt (fname, sol_lev[i], fmt='%.20lf')

rhs_lev = amg.get_rhs ()
for i in xrange(0, len (rhs_lev)):
    fname = "rhs." + repr(i)
    numpy.savetxt (fname, rhs_lev[i], fmt='%.20lf')

exit()

# output AMG setup data:
m_lev = amg.get_matrices ()
cf_lev = amg.get_cf_markers ()
s_lev = amg.get_s_markers ()
p_lev = amg.get_p_matrices ()

for i in xrange(0, len (m_lev)):
    fname = "A." + repr(i) + ".csr"
    t.ascii_write_to_csr_format (m_lev[i], fname, 0)
    
for i in xrange(0, len (p_lev)):
    fname = "P." + repr(i) + ".csr"
    t.ascii_write_to_csr_format (p_lev[i], fname, 0)
    
for i in xrange(0, len (s_lev)):
    fname = "S." + repr(i)
    numpy.savetxt (fname, s_lev[i], fmt='%d')

for i in xrange(0, len (cf_lev)):
    fname = "CF." + repr(i)
    numpy.savetxt (fname, cf_lev[i], fmt='%d')
