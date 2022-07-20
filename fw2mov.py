import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np
import os, sys, getopt, time

# progress bar - for looks and debugging
def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.3+
    count = len(it)
    def show(j):
        x = int(size*j/count)
        print("{}[{}{}] {}/{}".format(prefix, "#"*x, "."*(size-x), j, count), 
                end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("", flush=True, file=out)

# this was taken out of FUNWAVE-TVD plot_levee_solitary.py
def readETAData(num,numOfETA, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'eta_{0:05d}'.format(fileIndex)
    freeSurface = np.fromfile(fileName).reshape((nglob, mglob))
    return freeSurface

# this sure is the main function
def main(argv):
    # default vals for script
    dpi = 300   # dpi for figure
    wid = 18
    leg = 5     # curse len()... this is length
    dim = ''
    xmin = -1
    xmax = -1
    ymin = -1
    ymax = -1

    # parse *argv[] (is this securely done? no idea what's under the hood in c, so maybe...)
    logfilename = ''
    if len(argv) > 10:
        print('warning: new args ignored\n')
    try:
        opts, args = getopt.getopt(argv, "hi:o:d:w:l:",
            ["ifile=", "odir=", "dpi=", "wid=","len=", 
            "1d", "2d", "xmin=", "xmax=", "ymin", "ymax"])
    except getopt.GetoptError:
        print('fw2mov.py -i <inputlogfile.txt>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('fw2mov.py -i <inputlogfile.txt>')
            print('fw2mov should be in the working directory of the current model')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            if not ".txt" in arg:
                print("error: input log file not .txt")
                sys.exit(2)
            else:
                logfilename = arg
        elif opt in ("-d", "--dpi"):
            dpi = int(arg)
        elif opt in ("-w", "--wid"):
            wid = int(arg)
        elif opt in ("-l", "--len"):
            leg = int(arg)
        elif opt in ("--1d"):
            dim = "1d"
        elif opt in ("--2d"):
            dim = "2d"
        elif opt in ("--xmin"):
            xmin = int(arg)
        elif opt in ("--xmax"):
            xmin = int(arg)
        elif opt in ("--ymin"):
            ymin = int(arg)
        elif opt in ("--ymax"):
            ymax = int(arg)
    
    # parse log txt file
        # here's to hoping that the i/o of FUNWAVE-TVD doesn't get reworked
    print("--------------------")
    print("reading: " + logfilename)
    parse = ''
    with open(logfilename) as file:
        for line in file:
            if "input end" in line:
                break
            if not "---" in line:
                parse = [i for i in line.split() if i]
                for i in range(len(parse)):
                    if 'Mglob' in parse[i]:
                        Mglob = int(parse[i + 1])
                        print("parsed Mglob: " + str(Mglob))
                    elif "Nglob" in parse[i]:
                        Nglob = int(parse[i + 1])
                        print("parsed Nglob: " + str(Nglob))
                    elif "DX=" in parse[i]:
                        dx = float(parse[i + 1])
                        print("parsed dx: " + str(dx))
                    elif "DY=" in parse[i]:
                        dy = float(parse[i + 1])
                        print("parsed dy: " + str(dy))
                    elif "DEPTH_FILE" in parse[i]:
                        depthtype = "data"
                        depthfilepath = parse[i].split(':')[1]
                        print("parsed depth, depth type: " + depthtype)
                        print("parsed depthfilepath: "+ depthfilepath)
                    elif "DEPTH_FLAT" in parse[i]:
                        depthtype = "flat"
                        depthflat = float(parse[i+1])
                        print("parsed depth, depth type: flat (possibly slope)")
                        print("parsed depthflat: " + str(depthflat))
                    elif "SLP" in parse[i]:
                        depthtype = "slope"
                        slp = float(parse[i + 1])
                        print("parsed depth type: slope")
                        print("parsed slp: " + str(slp))
                    elif "Xslp" in parse[i]: 
                        xslp = float(parse[i + 1])
                        print("parsed xslp: " + str(xslp))
                    # elif "PLOT_INTV" in parse[i]:               #NOT SURE IF NEEDED...
                    #     plotint = float(parse[i + 1])
                    #     print("parsed plotint: " + str(plotint))
                    elif "RESULT_FOLDER" in parse[i]:
                        outputdir = parse[i].split(':')[1]
                        print("parsed output dir: "+ outputdir)
    print("closing: " + logfilename)
    print("--------------------")
    
    # dir work
    wdir = os.path.dirname(os.path.realpath(__file__))
    odir = os.path.join(wdir, outputdir)
    # print(odir) # debugging fragment
    
    # establish depth
    # TODO: pattern matching routine for depth types
    depth = np.full((Nglob, Mglob), -6) # TODO: temp
    
    # skim eta and sta
    print("skimming output folder")
    numeta = -1                     # correcting for ETA00000
    numsta = 0
    maxeta = 0
    for file in progressbar(os.listdir(os.fsencode(odir)), "", 40):
        filename = os.fsdecode(file)
        if "eta" in filename:
            numeta = numeta + 1
            fs = np.fromfile(os.path.join(odir, filename)).reshape((Nglob, Mglob))
            if (max(fs[1,]) >= maxeta):
                maxeta = max(fs[1,])
        elif "sta" in filename:
            numsta = numsta + 1
    print("found: {:d} eta files".format(numeta)) # why did i start c string formatting here? good question...
    print("found: {:d} sta files".format(numsta))
    print("found: {:.3f}m as max eta".format(maxeta))
    print("--------------------")
    
    print("producing figure") 
    fig = plt.figure(figsize=(wid,leg), dpi=dpi)
    ax = fig.add_subplot(1,1,1)
    
    if xmin < 0:
        xmin = 0
    if xmax < 0:
        xmax = Mglob
    if ymin < 0:
        ymin = 0
    if ymax < 0:
        ymax = Nglob

    writer = FFMpegWriter(fps = 10)

    # if 1d, 2d not specifed, default behavior
    if dim == '':
        if Nglob <= 3:
            dim = '1d'
        else:
            dim = '2d'

    if dim == '1d':
        with writer.saving(fig, 'mov.mp4', dpi):
            for i in progressbar(range(1,numeta + 1), "", 40):
                surf = readETAData(i, numeta, Mglob, Nglob, odir)
                x = np.asarray([float(xa) * dx for xa in range(Mglob)])
                c = int(Nglob / 2)

                plt.clf()
                plt.plot(x[xmin:xmax], surf[c,xmin:xmax], '-c', linewidth = 0.2)

                # Water Fill:
                plt.fill_between(x[xmin:xmax], depth[c,xmin:xmax], surf[c,xmin:xmax],
                                where = surf[c,xmin:xmax] > depth[c,xmin:xmax],
                                facecolor = 'cyan', interpolate = True)
                # Bottom Fill:
                plt.fill_between(x[xmin:xmax], min(depth[c,xmin:xmax])-1, depth[c,xmin:xmax], 
                                where= depth[c,xmin:xmax] > (depth[c,xmin:xmax]-2),       
                                facecolor = '0.35', hatch = 'X')

                plt.ylim((min(depth[c,:]) - 1, maxeta + 1))
                plt.xlim(xmin * dx, xmax * dx)
                writer.grab_frame()
            print("finishing...")
    else: 
        with writer.saving(fig, 'mov.mp4', dpi):
            print("need to implement")
    plt.close()

if __name__ == "__main__":
    startTime = time.time()
    main(sys.argv[1:])
    print("Runtime: {:f} seconds".format(time.time() - startTime))