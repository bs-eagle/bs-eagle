import bs
import os, sys, time, traceback, glob

class asserter (bs.bs_assert) :

    def __init__ (self, b = True, file = "", line = -1, cond_s = "") :
        bs.bs_assert.__init__(self) #, b, file, line)
        self.b = b
        self.file = file
        self.line = line
        self.cond_s = cond_s

    def set_var_list (self, var_list) :
        self.var_list = var_list

    def handle (self) :
        if (self.b == False) :
            print 'Assertion! Expression: %s \nFile: "%s", Line: %d \nValues: %s\n' % (self.cond_s, self.file, self.line, self.var_list)
            return True

        return True

class factory (bs.bs_assert_factory) :
    def __init__ (self) :
        bs.bs_assert_factory.__init__ (self)
        
    def make (self, b, file, line, cond_s) :
        return asserter (b, file, line, cond_s)

f = factory ()
a = bs.bs_assert ()
a.set_factory(f)

bos = bs.bs_bos_core
#bos.enable_fpu_exceptions ()

def get_param (index) :
    if (index >= len (sys.argv)) :
        return ""

    if (sys.argv[index] == "--verbose") :
        return get_param (index + 1)

    return sys.argv[index]
def get_named (name, default = False) :
    for arg in sys.argv :
        if (arg == name) :
            sys.argv.remove (arg)
            return True

    return default
def get_time () :
    return get_named ("--time")
def get_deep () :
    return get_named ("--deep")

is_deep = get_deep ()
is_time = get_time ()
is_verbose = get_named ("--verbose")
is_no_action = get_named ("--no")
is_print_files = get_named ("--print-files")
is_init_only = get_named ("--init-only", False)

def bs_amg () :
    amg                         = bs.bs_amg_solver.amg_solver_di ()
    amg_prop                    = bs.bs_amg_solver.amg_properties ()
    amg_prop.adaptive_threshold = 0
    amg_prop.update             = 0
    amg_prop.strength_threshold = 0.75

    amg.amg_prop                = amg_prop
    amg.coarser                 = bs.bs_amg_solver.coarse_pmis_2_di ()
    amg.interpolator            = bs.bs_amg_solver.interpolator_standart_2_di ()

    return amg

def one_phase_amg () :
    amg                         = bs.bs_amg_solver.amg_solver_di ()
    amg_prop                    = bs.bs_amg_solver.amg_properties ()
    amg_prop.adaptive_threshold = 0
    amg_prop.update             = 0
    amg_prop.strength_threshold = 0.75

    amg.amg_prop                = amg_prop
    amg.coarser                 = bs.bs_amg_solver.coarse_pmis_2_di ()
    amg.interpolator            = bs.bs_amg_solver.interpolator_standart_2_di ()

    return amg

def default_solver () :
    return bs.bs_base_linear_solvers.gmres_solver2_seq_di ()

def default_cfl_prec () :
    cpr     = bs.bs_cpr_prec.cpr_prec_seq_di ()
    csr     = bs.bs_bos_core.csr_ilu_cfl_prec_seq_di ()
    tsp     = bs.bs_bos_core.two_stage_prec_seq_di ()
    cpr.amg = bs_amg ()

    tsp.set_prec_1 (cpr)
    tsp.set_prec_2 (csr)

    return tsp

def default_prec () :
    cpr     = bs.bs_cpr_prec.cpr_prec_seq_di ()
    csr     = bs.bs_csr_ilu_prec.csr_ilu_prec_seq_di ()
    tsp     = bs.bs_bos_core.two_stage_prec_seq_di ()
    cpr.amg = bs_amg ()

    tsp.set_prec_1 (cpr)
    tsp.set_prec_2 (csr)

    return tsp

def ilu_prec () :
    csr = bs.bs_csr_ilu_prec.csr_ilu_prec_seq_di ()

    return csr

def ilu_cfl_prec () :
    csr = bs.bs_bos_core.csr_ilu_cfl_prec_seq_di ()

    return csr

def prepare_n_phase_solver (rs) :
    rs.jacobian.solver  = default_solver ()

    #if (rs.calc_model.params.get_by_name_i ("PREC_TYPE") == 1) :
    #    rs.jacobian.prec = ilu_prec ()
    #elif (rs.calc_model.params.get_by_name_i ("PREC_TYPE") == 2) :
    #    pass
    #else :
    #    rs.jacobian.prec = default_prec ()

    #if (rs.calc_model.params.get_by_name_b ("USE_CFL") == 1) :
    #    rs.jacobian.prec = ilu_cfl_prec ()
    #else :
    #    rs.jacobian.prec = ilu_prec ()
    if (rs.calc_model.params.get_by_name_b ("USE_CFL") == 1) :
        rs.jacobian.prec = default_cfl_prec ()
    else :
        rs.jacobian.prec = default_prec ()

def prepare_1_phase_solver (rs) :
    rs.jacobian.solver  = one_phase_amg ()
    if (rs.calc_model.params.get_by_name_i ("1_PHASE_PREC_TYPE") == 1) : 
        if (rs.calc_model.params.get_by_name_b ("USE_CFL") == 1) :
            rs.jacobian.prec = ilu_cfl_prec ()
        else :
            rs.jacobian.prec = ilu_prec ()
    else :
        rs.jacobian.prec = bs.bs_amg_solver.amg_solver_di ()

def simulate (input_name, is_no_action) :

    model_name = input_name.replace ("'", "")

    if (not is_no_action) :
        try :
            print ("Create reservoir simulator...")
            rs = bos.reservoir_simulator_di ()

            def on_post_read () :
                print ("on_post_read")

            rs.add_post_read_handler (on_post_read)

            print ("Load data from " + model_name)
            rs.init (model_name)

            print ("Set solver...")
            if (rs.calc_model.n_phases == 1) :
                prepare_one_phase_solver (rs)
            else :
                prepare_n_phase_solver (rs)

            if (not is_init_only) :
                print ("Simulate...")
                rs.simulate ()

        except:
            print ("Exception")
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            traceback.print_exception(exceptionType, exceptionValue, exceptionTraceback)
        print ("End")
    else :
        print ("No action for: " + model_name)

def process_file_list (file_list, is_time, is_no_action) :
    for file in file_list :
        if (is_time) :
            start = time.clock ()
            simulate (file, is_no_action)
            end = time.clock ()
            print ("Execution time: " + str (end - start))
        else :
            simulate (file, is_no_action)

def get_file_list (root, ext) :
    dir_list = glob.glob (os.path.join (root, "*"))
    file_list = []
    for dir in dir_list :
        if (os.path.isdir (dir)) :
            d = get_file_list (dir, ext)
            if (len (d) <> 0) :
                file_list.extend (d)

    d = glob.glob (os.path.join (root, ext))
    if (len (d) == 1) :
        file_list.extend ([d])

    return file_list

if (is_deep and is_print_files) :
    file_list = [file[0] for file in get_file_list ("./", "*.DATA")]
    print (file_list)

if (is_deep) :
    file_list = [file[0] for file in get_file_list ("./", "*.DATA")]
    process_file_list (file_list, is_time, is_no_action)

elif (len (sys.argv) == 2) :
    process_file_list ([sys.argv[1]], is_time, is_no_action)
else :
    if (is_verbose) :
        print (sys.argv)

    process_file_list (sys.argv[1:], is_time, is_no_action)

