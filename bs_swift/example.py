import bs
m = bs.mx.bcsr_matrix()
t = bs.mx.bcsr_matrix_tools ()
t.random_init(m,150,1,0.1,5)
amg = bs.swift.amg_solver()
amg.setup(m)