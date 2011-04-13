import bs
m = bs.mx.bcsr_matrix()
t = bs.mx.bcsr_matrix_tools ()
#fname=spe10.csr
#t.ascii_read_from_csr_format (m, fname)
t.random_init(m,50,1,0.1,5)
amg = bs.swift.amg_solver()
amg.setup(m)