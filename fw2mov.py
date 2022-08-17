import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np
import os, sys, getopt, time

# progress bar 
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

# readETAData and breadETAData are varients of code from FUNWAVE-TVD plot_levee_solitary.py
def readETAData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'eta_{0:05d}'.format(fileIndex)
    freeSurface = np.loadtxt(fileName)
    return freeSurface

# readETAData reading in binary
def breadETAData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'eta_{0:05d}'.format(fileIndex)
    freeSurface = np.fromfile(fileName).reshape((nglob, mglob))
    return freeSurface

# readMaskData reads... mask data
def readMaskData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'mask_{0:05d}'.format(fileIndex)
    freeSurface = np.loadtxt(fileName)
    return freeSurface

# readMaskData reading in binary
def breadMaskData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'mask_{0:05d}'.format(fileIndex)
    freeSurface = np.fromfile(fileName).reshape((nglob, mglob))
    return freeSurface

# readUData reads... udata
def readUData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'u_{0:05d}'.format(fileIndex)
    freeSurface = np.loadtxt(fileName)
    return freeSurface

# readUData reading in binary
def breadUData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'u_{0:05d}'.format(fileIndex)
    freeSurface = np.fromfile(fileName).reshape((nglob, mglob))
    return freeSurface

# readMaskData reads... vdata
def readVData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'v_{0:05d}'.format(fileIndex)
    freeSurface = np.loadtxt(fileName)
    return freeSurface

# readVData reading in binary
def breadVData(num, mglob, nglob, odir):
    fileIndex = num
    fileName = odir+'/'+'v_{0:05d}'.format(fileIndex)
    freeSurface = np.fromfile(fileName).reshape((nglob, mglob))
    return freeSurface

# wow! this sure is the main function
def main(argv):
    # default vals
    dpi = 300   # dpi for figure
    wid = 18    # fig wid
    leg = 5     # fig wid
    dim = ''    # 1d vs 2d, this gets updated towards the figure creation
    xmin = -1   # min/maxes are limits for figure, based on pure grid coor
    xmax = -1       # these get calculated before figure creation
    ymin = -1
    ymax = -1
    odir = ''
    outfiletype = 'ascii'   # output file i/o type, default to ascii (which is FUNWAVE default)
    fps = 15    # video fps
    mint = 0    # min time 
    maxt = sys.maxsize  # max time
    depthtype = ''
    ismask = False
    isquiver = False

    # parse *argv[] 
    try:
        opts, args = getopt.getopt(argv, "hi:o:",
            ["ifile=", "odir=", "dpi=", "wid=","len=", 
            "1d", "2d", "xmin=", "xmax=", "ymin=", "ymax=", 
            "fps=", "maxt=", "mint=", "arrow="])
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
        elif opt in ("--dpi"):
            dpi = int(arg)
        elif opt in ("--wid"):
            wid = int(arg)
        elif opt in ("--len"):
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
        elif opt in ("--fps"):
            fps = int(arg)
        elif opt in ("-maxt"):
            maxt = int(arg)
        elif opt in ("-mint"):
            mint = int(arg)
        elif opt in ("-o, --odir"):
            odir = arg
        elif opt in ("--arrow"):
            qint = int(arg)
            isquiver = True

    # parse log txt file
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
                        depthflat = float(parse[i+1])
                        if depthtype != "slope":
                            depthtype = "flat"
                            print("parsed depth, depth type: flat")
                        print("parsed depthflat: " + str(depthflat))
                    elif "BINARY" in parse[i]:
                        outfiletype = "binary"
                        print("parsed file i/o type: binary")
                    elif "SLP" in parse[i]:
                        depthtype = "slope"
                        slp = float(parse[i + 1])
                        print("parsed depth type: slope")
                        print("parsed slp: " + str(slp))
                    elif "Xslp" in parse[i]: 
                        xslp = float(parse[i + 1])
                        print("parsed xslp: " + str(xslp))
                    elif "PLOT_INTV" in parse[i]:               
                        plotint = float(parse[i + 1])
                        print("parsed plotint: " + str(plotint))
                    elif "RESULT_FOLDER" in parse[i]:
                        outputdir = parse[i].split(':')[1]
                        print("parsed output dir: "+ outputdir)
                    elif "OUT_MASK" in parse[i]:
                        if parse[i + 1] == 'T':
                            ismask = True
    print("closing: " + logfilename)
    print("--------------------")
    
    # dir work
    wdir = os.path.dirname(os.path.realpath(__file__))  # working directory
    if odir == '':      # check if output specified in argv[]
        odir = os.path.join(wdir, outputdir)
    else:
        print("warning: using output directory specifed by -o,--odir")
    
    # init depth arr
    depoutfile = os.path.join(odir, "dep.out")  # attempt to find dep.out
    if os.path.exists(depoutfile):  # if dep.out exists, use it to create depth data
        if outfiletype == "binary":
            depth = np.fromfile(depoutfile).reshape((Nglob, Mglob)) * -1
        else:
            depth = np.loadtxt(depoutfile) * -1
    elif depthtype == "data":       # else, resort to LOG.txt params
        depth = np.loadtxt(depthfilepath) * -1
    elif depthtype == "flat":
        depth = np.full((Nglob, Mglob), -depthflat)
    elif depthtype == "slope":
        depth = np.full((Nglob, Mglob), -depthflat)
        for i in range(int(xslp), Mglob + 1): # refactor this for accuracy
            depth[:,Mglob] = depth[:,Mglob] + (i - xslp) * slp

    # if 1d, 2d not specifed, default behavior
    if dim == '':
        if Nglob <= 3:
            dim = '1d'
        else:
            dim = '2d'
    
    # skim eta and sta
    print("skimming output folder")
    if os.path.exists(depoutfile):
        print("found: dep.out")     # print statement for clarity purposes
    numeta = 0                      # number of eta files
    numsta = 0                      # number of sta files
    nummask = 0                     # number of mask files
    maxeta = 0                      # max eta found
    mineta = 0                      # min eta found
    numu = 0
    numv = 0

    if outfiletype == 'binary': # adjusting file reading
        readfunc = lambda a, b, c, d : np.fromfile(os.path.join(a, b)).reshape((c, d))
    else:
        readfunc = lambda a, b, c, d : np.loadtxt(os.path.join(a, b))

    for file in progressbar(os.listdir(os.fsencode(odir)), "", 40): # read all files in output
        filename = os.fsdecode(file)
        if "eta_" in filename:
            numeta = numeta + 1
            fs = readfunc(odir, filename, Nglob, Mglob)
            if dim == '1d':
                if (max(fs[1,]) > maxeta):
                    maxeta = max(fs[1,])
                if (min(fs[1,]) < mineta):
                    mineta = min(fs[1,])
            else:
                for row in fs:
                    if (max(row) > maxeta):
                        maxeta = max(row)
                    if (min(row) < mineta):
                        mineta = min(row)
        elif "sta_" in filename:
            numsta = numsta + 1
        elif "mask_" in filename and not "mask9" in filename:
            nummask = nummask + 1
        elif "u_" in filename:
            numu = numu + 1
        elif "v_" in filename:
            numv = numv + 1
    print("found: {:d} eta files".format(numeta))
    print("found: {:d} mask files".format(nummask))
    #print("found: {:d} sta files".format(numsta))
    print("found: {:.3f}m as max eta".format(maxeta))
    print("found: {:.3f}m as min eta".format(mineta))
    if isquiver:
        print("found: {:d} u files".format(numu))
        if dim == '2d':
            print("found: {:d} v files".format(numv))

    # file checks:
    if not numeta == nummask:
        print("warning: masking is off - number of eta vs mask files not equal")
        ismask = False
    if isquiver:
        if not numeta == numu:
            print("warning: quiver is off - number of eta vs u files not equal")
            isquiver = False
        elif dim == '2d':
            if not numeta == numv:
                print("warning: quiver is off - number of eta vs v files not equal")
                isquiver = False
    print("--------------------")
    
    print("producing figure") 
    fig = plt.figure(figsize=(wid,leg), dpi=dpi)
    ax = fig.add_subplot(1,1,1)
    
    if xmin < 0:        # set up xmin, xmax, ymin, ymax
        xmin = 0
    if xmax < 0:
        xmax = Mglob
    if ymin < 0:
        ymin = 0
    if ymax < 0:
        ymax = Nglob

    writer = FFMpegWriter(fps = fps)

    if outfiletype == 'binary':
        readfunc = lambda a, b, c, d : breadETAData(a, b, c, d)
        maskfunc = lambda a, b, c, d : breadMaskData(a, b, c, d)
        ufunc = lambda a, b, c, d : breadUData(a, b, c, d)
        vfunc = lambda a, b, c, d : breadVData(a, b, c, d)
    else:
        readfunc = lambda a, b, c, d : readETAData(a, b, c, d)
        maskfunc = lambda a, b, c, d : readMaskData(a, b, c, d)
        ufunc = lambda a, b, c, d : readUData(a, b, c, d)
        vfunc = lambda a, b, c, d : readVData(a, b, c, d)

    if dim == '1d':
        with writer.saving(fig, 'mov_tint_{:d}_{:d}_xint_{:.0f}_{:.0f}.mp4'.format(mint, min(numeta,maxt), xmin, xmax), dpi):
            c = int(Nglob / 2)
            x = np.asarray([float(xa) * dx for xa in range(Mglob)])
            if isquiver:
                xu = np.asarray([float(xa) * dx * qint for xa in range(int(Mglob / qint))])
                yu = np.zeros(xu.size)
            for i in progressbar(range(mint, min(numeta, maxt)), "", 40):
                surf = readfunc(i, Mglob, Nglob, odir)
                if ismask:
                    mask = maskfunc(i, Mglob, Nglob, odir)
                    surf = np.ma.masked_where(mask == 0, surf)

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

                if isquiver:
                    ru = ufunc(i, Mglob, Nglob, odir)
                    u = np.zeros(xu.size)
                    for i in range(xu.size):
                        u[i] = ru[c,i * qint]
                    plt.quiver(xu[int(xmin / qint):int(xmax / qint)],
                        yu[int(xmin / qint):int(xmax / qint)], 
                        yu[int(xmin / qint):int(xmax / qint)],
                        u[int(xmin / qint):int(xmax / qint)],
                        scale = 10,
                        alpha = 0.4)

                plt.xlabel('Length (m)', fontsize = 12, fontweight = 'bold')
                plt.ylabel('Height (m)', fontsize = 12, fontweight = 'bold')
                plt.title("Time {:.5f} sec".format(plotint * i), fontsize = 12, fontweight = 'bold')

                writer.grab_frame()
            print("finishing...")
    else: 
        with writer.saving(fig, 'mov_tint_{:d}_{:d}_xint_{:.0f}_{:.0f}.mp4'.format(mint, min(numeta, maxt), xmin, xmax), dpi):
            x = np.asarray([float(xa)*dx for xa in range(Mglob)])
            y = np.asarray([float(ya)*dy for ya in range(Nglob)])
            if isquiver:
                xu = np.asarray([float(xa) * dx * qint for xa in range(int(Mglob / qint))])
                yu = np.asarray([float(ya) * dx * qint for ya in range(int(Nglob / qint))])
            for i in progressbar(range(mint, min(numeta, maxt)), "", 40):
                surf = readfunc(i, Mglob, Nglob, odir)
                if ismask:
                    mask = maskfunc(i, Mglob, Nglob, odir)
                    surf = np.ma.masked_where(mask == 0, surf)
                
                plt.clf()
                plt.pcolor(x[xmin:xmax], y[ymin:ymax], surf, cmap = "coolwarm")
                plt.xlabel('Y (m)', fontsize = 12, fontweight = 'bold')
                plt.ylabel('X (m)', fontsize = 12, fontweight = 'bold')
                plt.title("Time {:.5f} sec".format(plotint * i), fontsize = 12, fontweight = 'bold')
                plt.colorbar().set_label(r'$\eta$'+' (m)', rotation=90)
                plt.clim(mineta, maxeta)
                if isquiver:
                    ru = ufunc(i, Mglob, Nglob, odir)
                    rv = vfunc(i, Mglob, Nglob, odir)
                    u = np.zeros(ru.shape)
                    v = np.zeros(rv.shape)
                    np.copyto(u, ru)
                    np.copyto(v, rv)
                    plt.quiver(xu[int(xmin / qint):int(xmax / qint)],
                        yu[int(ymin / qint):int(ymax / qint)], 
                        u[int(xmin / qint):int(xmax / qint)],
                        v[int(ymin / qint):int(ymax / qint)],
                        scale = 10,
                        alpha = 0.2)
                writer.grab_frame()
            print("finishing...")
    plt.close()

if __name__ == "__main__":
    startTime = time.time()
    main(sys.argv[1:])
    print("runtime: {:f} seconds".format(time.time() - startTime))