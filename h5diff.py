import h5py, os, sys, unicodedata, glob, traceback
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

class make_diff :
    def __init__ (self, time, dates, d_params) :
        self.time = time
        self.dates = dates
        self.d_params = d_params

    def diff (self, index) :

        params = self.d_params

        s = "%d, %d" % (len (params[0][index]), len (self.dates[0]))
        assert len (params[0][index]) == len (self.dates[0]), s

        x = dict ()
        for i in xrange (0, len (params[0][index])) :
            p0 = params[0][index][i]
            p1 = params[1][index][i]
            date = self.dates[0][i]
            x[date] = [((p0 - p1) / max (p0, p1) if max (p0, p1) <> 0 else 0), p0, p1]

        for time in self.time :
            if time not in x :
                x[time] = [0, 0, 0]

        return x

def draw_diff (dir_name, well_name, data, name, range) :
    for i in range :
        title (dir_name + "/" + well_name + " (" + name[i] + ", diff)")
        d = [data[i][j][0] for j in xrange (0, len (data[i]))]
        plot (xrange (0, len (d)), d)
        savefig ("diff." + dir_name.replace ("./", "") + "_" + well_name + "." + name[i] + ".png")
        clf ()

def draw_data (dir_name, well_name, time, data, name, range) :
    for i in range :
        title (dir_name + "/" + well_name + " (" + name[i] + ", data)")
        lhs = [data[i][t][1] for t in time]
        rhs = [data[i][t][2] for t in time]
        plot (xrange (0, len (lhs)), lhs, xrange (0, len (rhs)), rhs)
        legend (('bs', 'v22'))
        savefig (dir_name + "_" + well_name + "." + name[i] + ".png")
        clf ()

def draw_acc (dir_name, well_name, time, data, name, range) :
    for i in range :
        title (dir_name + "/" + well_name + " (" + name[i] + ", acc)")
        lhs = [data[i][t][1] for t in time]
        rhs = [data[i][t][2] for t in time]
        plot (xrange (0, len (lhs)), lhs, xrange (0, len (rhs)), rhs)
        legend (('bs', 'v22'))
        savefig (dir_name + "_" + well_name + "." + name[i] + ".png")
        clf ()

def stat_diff (well_name, data, name, range, stat) :
    s = "Well name: " + well_name
    for i in range :

        diff = []
        line = []
        if (len (data[i]) == 0) :
            s += "\n\t%6s: no data" % (name[i])
            continue

        for j in xrange (0, len (data[i])) :
            t = data[i][j]
            if (len (diff) == 0) :
                diff = [t[0], t[0], t[0]]
                line = [0, j, j]
            else :
                diff[0] += t[0]
                if (fabs (t[0]) > fabs (diff[1])) :
                    diff[1] = t[0]
                    line[1] = j
                elif (fabs (t[0]) < fabs (diff[2])) :
                    diff[2] = t[0]
                    line[2] = j

        avg_diff = diff[0] / len (data[i])
        avg_diff_2 = avg_diff
        if (len (data[i]) - 1 > 0) :
            (diff[0] - diff[1]) / (len (data[i]) - 1)

        s += "\n\t%9s: avg: %9.6f (avg-wo-max: %9.6f), max: %9.6f (%04d), min: %9.6f (%04d)" % (name[i], avg_diff, avg_diff_2, diff[1], line[1], diff[2], line[2])

    stat.append (s)

def sum_total (time, data, range, total_data) :
    for i in range :
        try :
            if (i not in total_data) :
                total_data[i] = dict ()
                for t in time :
                    total_data[i][t] = [0, 0, 0]

            for t in time :
                total_data[i][t][1] += data[i][t][1]
                total_data[i][t][2] += data[i][t][2]
        except:
            print (data[i])
            raise

def draw_total (dir_name, time, data, name, range) :
    for i in range :
        title (dir_name + " (" + name[i] + ", total)")
        lhs = [data[i][t][1] for t in time]
        rhs = [data[i][t][2] for t in time]
        plot (xrange (0, len (lhs)), lhs, xrange (0, len (rhs)), rhs)
        legend (('bs', 'v22'))
        savefig (dir_name + "." + name[i] + ".png")
        clf ()

        if (get_verbose ()) :
            print (dir_name, name[i])

def diff_results (dir, lhs, rhs, what, model) :
    file = [h5py.File (lhs), h5py.File (rhs)]
    data = [file[0]["results"][what]["values"], file[1]["results"][what]["values"]]

    img = np.zeros ((len (data[0]), len (data[0][0])))
    data_len = len (data[0])
    if (what in ["pressure", "swat", "sgas"]) :
        for i in xrange (0, len (data[0])) :
            for j in xrange (0, len (data[0][i])) :
                p0 = data[0][i, j]
                p1 = data[1][i, j]
                img[i, j] = (p0 - p1)
            print (i, data_len)
    else:
        for i in xrange (0, len (data[0])) :
            for j in xrange (0, len (data[0][i])) :
                p0 = data[0][i, j]
                p1 = data[1][i, j]
                img[i, j] = ((p0 - p1) / max (p0, p1) if max (p0, p1) <> 0 else 0)
            print (i, data_len)

    title ("%s (%s)" % (dir, what))
    imshow (img, aspect='equal', extent=(0, 10000, 0, 10000))
    #imshow (img)
    cb=colorbar()
    #show ()
    savefig ("%s.%s.png" % (dir, what))
    clf ()

def draw_results (dir, lhs, rhs, what, model) :
    file = [h5py.File (lhs), h5py.File (rhs)]
    data = [file[0]["results"][what]["values"], file[1]["results"][what]["values"]]

    data_len = len (data[0][0])
    for l in xrange (0, data_len) :
        img = [np.zeros ((100, 100)), np.zeros ((100, 100))]
        j = 0
        k = 0
        for i in xrange (0, len (data[0])) :
            img[0][j, k] = data[0][i][l]
            img[1][j, k] = data[1][i][l]

            k += 1
            if (k == 100) :
                k = 0
                j += 1

        title ("%s (%s, lhs, %d)" % (dir, what, l))
        imshow (img[0])
        cb = colorbar ()
        savefig ("%s.%s.%d.lhs.png" % (dir, what, l))
        clf ()

        title ("%s (%s, rhs, %d)" % (dir, what, l))
        imshow (img[1])
        cb = colorbar ()
        savefig ("%s.%s.%d.rhs.png" % (dir, what, l))
        clf ()

        print (l, data_len)


def diff_wells (dir, lhs, rhs, action, what, model) :

    if (action == "diff" and skip_model (what, dir)) :
        return
    elif (action == "draw" and skip_model (model, dir)) :
        return

    file = [h5py.File (lhs), h5py.File (rhs)]
    wells = [file[0]["wells"], file[1]["wells"]]

    assert (len (wells[0]) == len (wells[1]))

    name = ["oil", "water", "gas", "liquid", "inj_oil", "inj_water", "inj_gas",
            "", "", "", "", "", "", "", 
            "", "",
            "acc_oil", "acc_water", "acc_gas", "acc_liquid", "acc_inj_oil", "acc_inj_water", "acc_inj_gas"]

    time = file[0]["results"]["pressure"]["dates"][0]

    stat = []
    total_data = dict ()
    for well in wells[0] :

        well_name = str (well).decode ("ascii", "ignore").encode ("ascii")

        d_params = [wells[0][well]["d_params"], wells[1][well]["d_params"]]
        dates = [wells[0][well]["dates"], wells[1][well]["dates"]]

        for i in xrange (0, len (dates[0])) :
            if (dates[0][i] <> dates[1][i]) :
                print ("DATES mismatch", dates[0][i], dates[1][i], well_name)
                #assert (dates[0][i] == dates[1][i])

        md = make_diff (time, dates, d_params)

        data = [md.diff (0), md.diff (1), md.diff (2), md.diff (3),
                md.diff (4), md.diff (5), md.diff (6),
                [], [], [], [], [], [], [],
                [], [],
                md.diff (16), md.diff (17), md.diff (18), md.diff (19), 
                md.diff (20), md.diff (21), md.diff (22)]

        if (action == "draw") :
            if (what == "diff") :
                draw_diff (dir, well_name, data, name, xrange (16, 23))
            elif (what == "data") :
                draw_data (dir, well_name, time, data, name, xrange (0, 4))
            elif (what == "acc") :
                draw_acc (dir, well_name, time, data, name, xrange (16, 23))
            elif (what == "total") :
                sum_total (time, data, xrange (16, 23), total_data)
                sum_total (time, data, xrange (0, 7), total_data)
            else :
                #draw_diff (dir, well_name, data, name, data_count)
                #draw_data (dir, well_name, data, name, data_count)
                pass
        elif (action == "diff") :
            stat_diff (well_name, data, name, xrange (16, 23), stat)

    if (action == "diff") :
        print ("Dir: %s" % (dir))
        for s in stat :
            print (s)
        print ("-------------------")
    elif (action == "draw" and what == "total") :
        draw_total (dir, time, total_data, name, xrange (16, 23))
        draw_total (dir, time, total_data, name, xrange (0, 7))

def get_param (index) :
    if (index >= len (sys.argv)) :
        return ""

    if (sys.argv[index] == "--verbose") :
        return get_param (index + 1)

    return sys.argv[index]

def get_verbose () :
    for arg in sys.argv :
        if (arg == "--verbose") :
            return True

    return False

def skip_model (model, dir) :
    if (get_verbose ()) :
        print (model, dir)

    if (len (model) <> 0 and model != dir + "/") :
        return True 

    return False

def get_file_list (root) :
    dir_list = glob.glob (os.path.join (root, "*"))
    file_list = []
    for dir in dir_list :
        if (os.path.isdir (dir)) :
            d = get_file_list (dir)
            if (len (d) <> 0) :
                file_list.extend (d)

    d = glob.glob (os.path.join (root, "*.h5"))
    if (len (d) == 2) :
        file_list.extend ([d])

    return file_list

action  = get_param (1)
what    = get_param (2)
model   = get_param (3)
verbose = get_verbose ()

file_list = get_file_list ("./")

for file in file_list :
    try:
        dir = os.path.dirname (file[0])
        y = [file[1], file[0]]
        if (file[0].find ("results.h5") != -1) :
            y = [file[0], file[1]]

        if (action == "diff") :
            if (what in ["pressure", "swat", "sgas"] and not skip_model (model, dir)) :
                diff_results (dir, y[0], y[1], what, model)
            else :
                diff_wells (dir, y[0], y[1], action, what, model)
        elif (action == "draw") :
            if (what in ["swat", "sgas"] and not skip_model (model, dir)) :
                draw_results (dir, y[0], y[1], what, model)
            else :
                diff_wells (dir, y[0], y[1], action, what, model)
        else:
            diff_wells (dir, y[0], y[1], action, what, model)
    except:
        print ("Diff failed: " + dir)
        #traceback.print_exc ()
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        traceback.print_exception(exceptionType, exceptionValue, exceptionTraceback)
