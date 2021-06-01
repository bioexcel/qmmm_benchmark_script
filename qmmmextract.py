import sys
import os
import numpy as np
import argparse

# returns values from CP2K output
def findvalueforpattern(filename, pattern):
    with open(filename, "r") as f:
        for line in f:
            if pattern in line:
                length = len(line.split())
                return line.split()[length-1]
    return False

# get the time from CP2K output, function for a specific function, 
# or total for total run time
def get_time(filename, function):
    if function=='total':
        pattern = 'CP2K   '
        column = 6
    else:
        pattern = function
        column = 4
    time = 0.0
    with open(filename, "r") as f:
        for line in f:
            if pattern in line:
                if len(line.split())==7:
                    time = float(line.split()[column])
                    break
    if time==0.0:
        print("time value not found in ", filename)
    return time

# sort self times in order
def sortselftime(key):
    return float(key['stime'][0])


# find the n subroutines with top self time
def top_self_time(filename, n):
    selftimes=[]
    with open(filename, "r") as f:
        for line in f:
            if 'CP2K   ' in line:
                break
        for line in f:
            if '-------' in line:
                break
            item = {}
            item['name'] = line.split()[0]
            item['stime'] = [float(line.split()[4])]
            selftimes.append(item)

        selftimes.sort(key=sortselftime, reverse=True)
        del selftimes[int(n):]
    return selftimes
 
# return standard error
def std_err(data):
    data = np.asarray(data)
    if np.all((data==0)):
        return 0
    error = np.std(data, ddof=1) / np.sqrt(np.size(data))
    return error

       
# append matching self times together
def appendselftimes(selftimes1, selftimes2):
    for item1 in selftimes1:
        present = False
        for item2 in selftimes2:
            if item1['name'] == item2['name']:
                item1['stime'].append(item2['stime'][0])
                present = True
   
    return selftimes1

# average self times
def aveselftimes(selftime):
    for item in selftime:
        average = sum(item['stime'])/len(item['stime'])
        #error = std_err(item['stime'])
        item['average'] = average
        #item['error'] = error
    return selftime

# sorting functions for the read in file data
def sortcores(key):
    return int(key['cores'])
def sortthreads(key):
    return int(key['threads'])
def sortsteps(key):
    return int(key['steps'])
def sortprojects(key):
    return key['project']

# format output
def fmt(x):
    return format(x, '0.5f')

# extract data for CP2K output files ending in .log
# report top n most costly subroutines
def extract_data(function, n):
    path = os.getcwd()
    fname = []
    for root,d_names,f_names in os.walk(path):
        for f in f_names:
            if not os.path.basename(root)=='restart':
                fname.append(os.path.join(root, f))

    filedata = []
    for file in fname:
        if file.endswith(".log"):
            print("Found CP2K logfile: " + file)
            fileinfo = {}
            steps = findvalueforpattern(file, 'Number of time steps')
            threads = findvalueforpattern(file, 'Number of threads for this process')
            processes = findvalueforpattern(file, 'Total number of message passing processes')
            project = findvalueforpattern(file, 'GLOBAL| Project name')
            cores = int(threads) * int(processes)
            time = get_time(file, function)
            # if steps not in output set to 1
            if not findvalueforpattern(file, 'Number of time steps'):
                fileinfo['steps'] = "1"
            else:
                fileinfo['steps'] = steps
            fileinfo['threads'] = threads
            fileinfo['project'] = project
            fileinfo['cores'] = str(cores)
            fileinfo['time'] = [time]
            fileinfo['selftimes'] = top_self_time(file, n)

            # check if same run already exists
            found = False
            for cp2kout in filedata:
                # if exists we append the times found for aveaging later
                if (cp2kout["steps"] == steps and cp2kout["threads"]==threads and cp2kout['cores']==str(cores)) and cp2kout['project']==project:
                    cp2kout["time"].append(time)
                    cp2kout["selftimes"] = appendselftimes(cp2kout['selftimes'], top_self_time(file, n))
                    found = True
            # if not exists we add a new entry
            if found==False:
                filedata.append(fileinfo)
    return filedata


def extract(function="total", n=5):

    # extract the data from files
    filedata = extract_data(function, n)

    # find all unique values of cores, threads, steps, projects
    corevals = set()
    threadvals = set()
    stepvals = set()
    projects = set()
    for cp2kout in filedata:
        threadvals.add(cp2kout['threads'])
        stepvals.add(cp2kout['steps'])
        projects.add(cp2kout['project'])
    # set up output files for the data   
    for threads in threadvals:
        for steps in stepvals:
            for proj in projects:
                output = proj + "_" + function + "_" + steps + "steps_" + threads + "threads.out"
                barfile = proj + "_topcalls_" + steps + "steps_" + threads + "threads.out"
                # remove files if they exist
                if os.path.exists(output):
                    os.remove(output)
                if os.path.exists(barfile):
                    os.remove(barfile)
                uniquesubroutinenames = set()
                # find unique top subroutines
                for cp2kout in filedata:
                    if cp2kout['steps'] == steps and cp2kout['threads']==threads and cp2kout['project']==proj:
                        for values in cp2kout['selftimes']:
                            uniquesubroutinenames.add(values['name'])
                if len(uniquesubroutinenames) !=0:
                    bar = open(barfile, "w")
                    bar.write("cores\t")
                    for names in uniquesubroutinenames:
                        bar.write(names + "\t")
                    bar.close()

    # sort data
    filedata.sort(key=sortcores)
    filedata.sort(key=sortthreads)
    filedata.sort(key=sortsteps)
    filedata.sort(key=sortprojects)

    #write out time info
    for cp2kout in filedata:
        output = cp2kout['project'] + "_" + function + "_" + cp2kout['steps'] + "steps_" + cp2kout['threads'] + "threads.out"
        if os.path.exists(output):
            of = open(output, "a")
        else:
            of = open(output, "w")
            of.write("cores     time     error\n")
        of.write(cp2kout['cores'] + "\t" + fmt(sum(cp2kout['time'])/len(cp2kout['time'])) + "\t" + fmt(std_err(cp2kout['time'])) + "\n")
        of.close()

    # write out subroutine top times
    for cp2kout in filedata:
         # average the selftimes 
        cp2kout['selftimes'] = aveselftimes(cp2kout['selftimes'])
        totaltimeaverage = sum(cp2kout["time"])/len(cp2kout["time"])
        barfile = cp2kout['project'] + "_topcalls_" + cp2kout['steps'] + "steps_" + cp2kout['threads'] + "threads.out"
        bar = open(barfile, "r")
        # read first line to find subroutnie list
        bar.seek(0)
        subroutinenames = bar.readline()
        subroutinenames = subroutinenames.split()
        del subroutinenames[0]
        bar.close()
        bar = open(barfile, "a")
        bar.write("\n" + cp2kout['cores'] + "\t")
        for names in subroutinenames:
            exists=False
            for value in cp2kout['selftimes']:
                if value['name'] == names:
                    bar.write(fmt(value['average']/totaltimeaverage) + "\t")
                    exists=True
            if exists==False:
                bar.write("0.00   \t")
        bar.write("\n")


def plotdata(filename, plottype="time"):

    if not os.path.exists(filename) or filename=='all':
        sys.exit("filename doesn't exist")

    if (plottype=='time'):
        plotscript = "time.plt"
        ylabel = "Time [s]"
        xlabel = "Cores"
    elif (plottype=='speedup'):
        plotscript = "speedup.plt"
        ylabel = "Speed up"
        xlabel = "Number of nodes"
    elif (plottype=='bar'):
        plotscript = "bar.plt"
        ylabel = "Fraction of total time"
        xlabel = "Cores"

    f = open(plotscript, "w")
    f.write('#plot script for ' + plottype + ' \n')
    f.write('set key spacing 1.5 \n')
    f.write('set key top right \n')

    if plottype=='bar':
        f.write('set style data histograms \n')
        f.write('set style histogram rowstacked \n')
        f.write('set boxwidth 0.5 relative \n')
        f.write('set style fill solid 1.0  \n')
        f.write('set key at 9,16 \n')
        f.write('set key spacing 1.1 \n')
        f.write('set key font "Helvetica, 15" \n')
        f.write('set rmargin 40 \n')
    f.write('set xlabel font "Helvetica, 26" \n')
    f.write('set ylabel font "Helvetica, 26" offset -2 \n')
    f.write('set ytics font "Helvetica, 20" \n')
    f.write('set xtics font "Helvetica, 20" \n')
    f.write('set linetype  1 lc rgb "dark-violet" lw 4 \n')
    f.write('set linetype  2 lc rgb "#0011a7" lw 4 \n')
    f.write('set linetype  3 lc rgb "#ff0101" lw 4 \n')
    f.write('set linetype  6 lc rgb "#9ed323" lw 4 \n')
    f.write('set linetype  7 lc rgb "#fea655" lw 4 \n')
    f.write('set linetype  5 lc rgb "#70cce1" lw 4 \n')
    f.write('set linetype  4 lc rgb "#0b8900" lw 4 \n')
    f.write('set linetype  8 lc rgb "black"   lw 4 \n')
    f.write('set linetype  9 lc rgb "gray50"  lw 4 \n')
    f.write('set size ratio 0.5 \n')
    f.write('set lmargin 12 \n')
    f.write('set xlabel "' + xlabel + '"' + '\n')
    f.write('f(x) = x \n')
    f.write('set ylabel "' + ylabel +'"' + '\n')
    f.write('set terminal postscript enhanced color eps background rgb "white" \n')
    f.write('set output sprintf("%s-' + plottype + '.eps", filename) \n')
    if plottype=='time':
        f.write('unset key \n')
        f.write('plot for [i=2:*:2] filename u 1:i w lp lt i title columnheader, for [i=2:*:2] filename u 1:i:i+1 w yerrorbars notitle ps 2 pt 2 lt i/2+1 lw 4 \n')
    if plottype=='speedup':
        f.write('unset key \n')
        f.write("""cpn = system("awk 'FNR == 2 {print $1}' """ + filename +'") \n')
        f.write("""first = system("awk 'FNR == 2 {print $2}' """ + filename +'") \n')
        f.write('plot filename u ($1/cpn):((first/$2)) w lp lw 4 ps 2 pt 2 lt 2, f(x) lw 2 lt 9 \n')
    if plottype=='bar':
        f.write('set key noenhanced \n')
        f.write('plot for [i=2:*:1] filename u i:xticlabels(1) lt i title columnheader \n')
    f.write('exit \n')

    f.close()
    path = os.getcwd()

    if filename=='all':
        for file in os.listdir(path):
            if file.endswith(".out"):
                print(os.path.join(path, file))
                gnuplotcommand = """gnuplot -e "filename='""""" + path + file + "'" + '"' + " " + plotscript
                print(gnuplotcommand)
                os.system(gnuplotcommand)
    
    gnuplotcommand = """gnuplot -e "filename='""""" + filename + "'" + '"' + " " + plotscript
    os.system(gnuplotcommand)


def readdata(filename):
    data=[]
    cores=[]
    error=[]
    with open(filename, "r") as f:
        f.readline()
        for line in f:
            cores.append(line.split()[0])
            data.append(float(line.split()[1]))
            error.append(float(line.split()[2]))
    return cores,data,error



#def compareplot(file1, file2, cores):
    # paste files together 
    # sensible column names
    # plot data

def gettimeperstep(name, function, val1, val2, threads):



    datafilename1 = name + "_" + function + "_" + str(val1) + "steps_" + str(threads) + "threads.out"
    datafilename2 = name + "_" + function + "_" + str(val2) + "steps_" + str(threads) + "threads.out"
    difffilename = name + "_" + function + "_" + str(val1) + "-"+ str(val2) + "steps_" + str(threads) + "threads.out"

    if not os.path.isfile(datafilename1):
        sys.exit("filename " + datafilename1 + " doesn't exist")
    if not os.path.isfile(datafilename2):
        sys.exit("filename " + datafilename2 + " doesn't exist")

    of = open(difffilename, "w")
    of.write("cores     time     error\n")
    difference = int(val1) - int(val2)
    cores1,data1,error1 = readdata(datafilename1)
    cores2,data2,error2 = readdata(datafilename2)
    if cores1 == cores2:
        print("data matches")
    for i in range(0,len(data1)):
        of.write(cores1[i] + "\t" + fmt((data1[i] - data2[i])/difference) + "\n")


if __name__=='__main__':

    my_parser = argparse.ArgumentParser(prog="qmmmextract.py",
                                        usage="%(prog)s mode inputname option", 
                                        description="")
    my_parser.add_argument("mode", metavar="mode", type=str, 
                           help="The mode of operation")

    my_parser.add_argument("--filename", metavar="filename", type=str, 
                           help="The input option")
    my_parser.add_argument("--function", metavar="function", type=str, 
                           help="The input option")
    my_parser.add_argument("--notop", metavar="notop", type=str, 
                           help="The input option")

    my_parser.add_argument("--plottype", metavar="opt", type=str, 
                           help="The extra option")

    my_parser.add_argument("--name", metavar="notop", type=str, 
                           help="The input option")
    my_parser.add_argument("--maxsteps", metavar="notop", type=str, 
                           help="The input option")
    my_parser.add_argument("--minsteps", metavar="notop", type=str, 
                           help="The input option")
    my_parser.add_argument("--threads", metavar="notop", type=str, 
                           help="The input option")
 
    args = my_parser.parse_args()
    mode = args.mode
    
    if mode=='extract':

        if args.function == None:
            print("No --function specified, doing total run time as default")
            function='total'
        else:
            function = args.function
        if args.notop == None:
            print("No --notop specified, doing the top 5 most costly subroutines")
            n = 5
        else:
            n = args.notop
            print("Finding the " + n + " most costly subroutines")
        extract(function, n)
    elif mode=='plot':


        if args.filename == None:
            sys.exit("Filename (--filename) required for plot option")
        else:
            filename = args.filename
        if args.plottype == None:
            print("No --plottype specified, doing time plot as default plot")
        else:
            plottype = args.plottype
        plotdata(filename, plottype)
    elif mode=='timestep':
        
        if args.name == None:
            sys.exit("CP2K project name (--name) required for timestep option")
        else:
            name = args.name
        if args.function == None:
            sys.exit("Subroutine name (--function) required for timestep option")
        else:
            function = args.function

        if args.maxsteps == None:
            sys.exit("Max steps (--maxsteps) required for timestep option")
        else:
            maxsteps = args.maxsteps
        if args.minsteps == None:
            sys.exit("Min steps (--minsteps) required for timestep option")
        else:
            minsteps = args.minsteps
        if args.threads == None:
            sys.exit("Number of threads (--threads) required for timestep option")
        else:
            threads = args.threads

        threads = threads.split(',')
        for t in threads:
            gettimeperstep(name,function,maxsteps,minsteps,t)
    #elif mode==compare:

    else:
        sys.exit("mode not recognised")
