import sys
import os
import numpy as np
import argparse

def findvalueforpattern(filename, pattern):
    with open(filename, "r") as f:
        for line in f:
            if pattern in line:
                length = len(line.split())
                return line.split()[length-1]
    return False

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

def sorttime(key):
    return float(key['stime'][0])

def top_self_time(filename):
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

        selftimes.sort(key=sorttime, reverse=True)
        del selftimes[10:]
    return selftimes
 
def std_err(data):
    data = np.asarray(data)
    error = np.std(data, ddof=1) / np.sqrt(np.size(data))
    return error

       

def appendselftimes(selftimes1, selftimes2):
    for item1 in selftimes1:
        present = False
        for item2 in selftimes2:
            if item1['name'] == item2['name']:
                item1['stime'] = item1['stime']
                item1['stime'].append(item2['stime'][0])
                present = True
   
    return selftimes1

def aveselftimes(selftime):
    for item in selftime:
        average = sum(item['stime'])/len(item['stime'])
        error = std_err(item['stime'])
        item['average'] = average
        item['error'] = error
    return selftime

def getuniquenames(filedata):
    uniquenames = set()
    for item in filedata:
        for values in item['selftimes']:
            uniquenames.add(values['name'])
    return uniquenames

def sortcores(key):
    return int(key['cores'])
def sortthreads(key):
    return int(key['threads'])
def sortsteps(key):
    return int(key['steps'])
def sortprojects(key):
    return key['project']


def fmt(x):
    return format(x, '0.5f')

def extract_data(function):
    path = os.getcwd()
    fname = []
    for root,d_names,f_names in os.walk(path):
        for f in f_names:
            if not os.path.basename(root)=='restart':
                fname.append(os.path.join(root, f))

    filedata = []
    for file in fname:
        if file.endswith(".log"):
            print(file)
            fileinfo = {}
            steps = findvalueforpattern(file, 'Number of time steps')
            threads = findvalueforpattern(file, 'Number of threads for this process')
            processes = findvalueforpattern(file, 'Total number of message passing processes')
            project = findvalueforpattern(file, 'GLOBAL| Project name')
            cores = int(threads) * int(processes)
            time = get_time(file, function)
            if not findvalueforpattern(file, 'Number of time steps'):
                fileinfo['steps'] = "1"
            else:
                fileinfo['steps'] = steps
            fileinfo['threads'] = threads
            fileinfo['project'] = project
            fileinfo['cores'] = str(cores)
            fileinfo['time'] = [time]
            fileinfo['selftimes'] = top_self_time(file)
            found = False
            for item in filedata:
                if (item["steps"] == steps and item["threads"]==threads and item['cores']==str(cores)) and item['project']==project:
                    item["time"].append(time)
                    item["selftimes"] = appendselftimes(item['selftimes'], top_self_time(file))
                    found = True
            if found==False:
                filedata.append(fileinfo)
    return filedata


def extract(function="total"):

    filedata = extract_data(function)

    uniquenames = getuniquenames(filedata)
    corevals = set()
    threadvals = set()
    stepvals = set()
    projects = set()
    for item in filedata:
        threadvals.add(item['threads'])
        stepvals.add(item['steps'])
        projects.add(item['project'])
        item['selftimes'] = aveselftimes(item['selftimes'])
    for threads in threadvals:
        for steps in stepvals:
            for proj in projects:
                output = proj + "_" + function + "_" + steps + "steps_" + threads + "threads.out"
                barfile = proj + "_topcalls_" + steps + "steps_" + threads + "threads.out"
                if os.path.exists(output):
                    os.remove(output)
                of = open(output, "w")
                bar = open(barfile, "w")
                for names in uniquenames:
                    bar.write(names + "\t")
                of.write("cores     time     error\n")
                of.close()

    filedata.sort(key=sortcores)
    filedata.sort(key=sortthreads)
    filedata.sort(key=sortsteps)
    filedata.sort(key=sortprojects)

    #write out time info
    for item in filedata:
        output = item['project'] + "_" + function + "_" + item['steps'] + "steps_" + item['threads'] + "threads.out"
        of = open(output, "a")
        of.write(item['cores'] + "\t" + fmt(sum(item['time'])/len(item['time'])) + "\t" + fmt(std_err(item['time'])) + "\n")

    # write out subroutine top times
    for item in filedata:
        barfile = item['project'] + "_topcalls_" + item['steps'] + "steps_" + item['threads'] + "threads.out"
        bar = open(barfile, "a")
        bar.write("\n" + item['cores'] + "\t")
        for names in uniquenames:
            exists=False
            for value in item['selftimes']:
                if value['name'] == names:
                    bar.write(str(value['average']) + "\t")
                    exists=True
            if exists==False:
                bar.write("?   \t")
        bar.write("\n")


def plotdata(filename, plottype="time"):

    if not os.path.exists(filename):
        sys.exit("filename doesn't exist")

    if (plottype=='time'):
        plotscript = "time.plt"
        ylabel = "Time [s]"
    if (plottype=='speedup'):
        plotscript = "speedup.plt"
        ylabel = "Speed up"

    f = open(plotscript, "w")
    f.write('#plot script for speed up \n')
    f.write('set key spacing 1.5 \n')
    f.write('set key top right \n')
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
    f.write('set xlabel "Cores" \n')
    f.write('f(x) = x \n')
    f.write('set ylabel "' + ylabel +'"' + '\n')
    f.write('set terminal postscript enhanced color eps \n')
    f.write('set output sprintf("%s.eps", filename) \n')
    if plottype=='time':
        f.write('plot for [i=2:*:2] filename u 1:i w lp lt i title columnheader, for [i=2:*:2] filename u 1:i:i+1 w yerrorbars notitle ps 2 pt 2 lt i/2+1 lw 4 \n')
    if plottype=='speedup':
        f.write('first(x) = ($0 > 0 ? base : base = x) \n')
        f.write('plot for [i=2:*] filename u 1:(first($i),base/$i) w lp lw 4 ps 2 pt 2 \n')
    f.write('exit \n')

    f.close()
    path = os.getcwd()

    '''
    for file in os.listdir(path):
        if file.endswith(".out"):
            print(os.path.join(path, file))
            gnuplotcommand = """gnuplot -e "filename='""""" + path + file + "'" + '"' + " " + plotscript
            print(gnuplotcommand)
            os.system(gnuplotcommand)
    '''
    gnuplotcommand = """gnuplot -e "filename='""""" + filename + "'" + '"' + " " + plotscript
    print(gnuplotcommand)
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


def gettimeperstep(name, function, val1, val2, threads):

    datafilename1 = name + "_" + function + "_" + str(val1) + "steps_" + str(threads) + "threads.out"
    datafilename2 = name + "_" + function + "_" + str(val2) + "steps_" + str(threads) + "threads.out"
    difffilename = name + "_" + function + "_" + str(val1) + "-"+ str(val2) + "steps_" + str(threads) + "threads.out"

    of = open(difffilename, "w")
    of.write("cores     time     error\n")
    difference = val1 - val2
    cores1,data1,error1 = readdata(datafilename1)
    cores2,data2,error2 = readdata(datafilename2)
    print (data1)
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

    my_parser.add_argument("inputname", metavar="in", type=str, 
                           help="The input option")

    my_parser.add_argument("--plottype", metavar="opt", type=str, 
                           help="The extra option")


    args = my_parser.parse_args()
    mode = args.mode
    
    if mode=='extract':
        function = args.inputname
        extract(function)
    elif mode=='plot':
        filename = args.inputname
        plottype = args.plottype
        plotdata(filename, plottype)
    elif mode==tperstep:
        gettimeperstep("MQAE","total",6,1,1)
    else:
        sys.exit("mode not recognised")
