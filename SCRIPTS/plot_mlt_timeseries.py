import numpy as np
import glob
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import re
import datetime as dt
import argparse
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

# def class runid
class run(object):
    def __init__(self, runid, cisf, sf):
        # parse dbfile
        self.runid, self.name, self.line, self.color = parse_dbfile(runid,cisf)
        self.sf = sf

    def load_time_series(self, cfile, cvar):
        # need to deal with mask, var and tag
        # need to do with the cdftools unit -> no unit !!!!
        # define time variable
        ctime = 'time_centered'

        # define unit
        nf = len(cfile)
        df=[None]*nf
        for kf,cf in enumerate(cfile):
            try:
                ncid    = nc.Dataset(cf)

                ncvtime = ncid.variables[ctime]
                if 'units' in ncvtime.ncattrs():
                    cunits = ncvtime.units
                else:
                    cunits = "seconds since 1900-01-01 00:00:00"
    
                # define calendar
                if 'calendar' in ncvtime.ncattrs():
                    ccalendar = ncvtime.calendar
                else:
                    ccalendar = "noleap"
                time = nc.num2date(ncid.variables[ctime][:].squeeze(), cunits, ccalendar)
        
                # convert to proper datetime object
                if isinstance(time,(list,np.ndarray)):
                    ntime = time.shape[0]
                else:
                    ntime = 1
    
                timeidx=[None]*ntime
                for itime in range(0, ntime):
                    if isinstance(time,(list,np.ndarray)):
                        timeidx[itime] = np.datetime64(time[itime],'us')
                    else:
                        timeidx[itime] = np.datetime64(time,'us')
        
                # build series
                cnam=get_varname(cf,cvar)
                df[kf] = pd.Series(ncid.variables[cnam][:].squeeze() * self.sf, index = timeidx, name = self.name)

            except Exception as e: 
                print 'issue in trying to load file : '+cf
                print e
                sys.exit(42) 


        # build dataframe
        self.ts   = pd.DataFrame(pd.concat(df)).sort_index()
        self.mean = self.ts[self.name].mean()
        self.std  = self.ts[self.name].std()
        self.min  = self.ts[self.name].min()
        self.max  = self.ts[self.name].max()

    def __str__(self):
        return 'runid = {}, name = {}, line = {}, color = {}'.format(self.runid, self.name, self.line, self.color)

def get_name(regex,varlst):
    revar = re.compile(r'\b%s\b'%regex,re.I)
    cvar  = filter(revar.match, varlst)
    if (len(cvar) > 1):
        print regex+' name list is longer than 1 or 0; error'
        print cvar[0]+' is selected'
    if (len(cvar) == 0):
        print 'no match between '+regex+' and :'
        print varlst
        sys.exit(42)
    return cvar[0]

def get_varname(cfile,cvar):
    ncid   = nc.Dataset(cfile)
    cnam=get_name(cvar,ncid.variables.keys())
    ncid.close()
    return cnam

#=============================== obs management =================================
class obs:
    def __init__(self, cobs, disf):
        # parse dbfile
        self.name= cobs
        self.mkr = self.load_obsmkr() # char
        self.dta = self.load_obsdta() # dic{isf,(dta err)}
        self.clr = disf               # dic{isf, clr}

    def load_obsmkr(self):
        try:
            lstyle=False
            with open('obs.sty') as fid:
                for cline in fid:
                    att=cline.split('|')
                    if att[0].strip() == self.name:
                        cpltmarker  = att[1].strip()
                        lstyle=True
            if not lstyle:
                print runid+' not found in obs.sty'
                raise Exception
    
        except Exception as e:
            print e
            sys.exit(42)
    
        return cpltmarker

    def load_obsdta(self):
        obs_isf={}
        with open(self.name+'.txt', 'r') as obsfile:
            for line in obsfile:
                if line.split()[0] == 'EOF':
                    break
                isfname    = line.split()[0]
                isfmlt_val = float(line.split()[1])
                isfmlt_err = float(line.split()[2])
                obs_isf[isfname]=(isfmlt_val, isfmlt_err)
        return obs_isf
    
def find_key(char, fid):
    for cline in fid:
        lmatch = re.findall(char, cline) 
        if (lmatch) :
            return cline.rstrip().strip('\n').split(' ')[-1]
#================================================================================

def load_isfclr():
    try:
        lstyle=False
        clr={}
        with open('isf.sty') as fid:
            for cline in fid:
                att=cline.split('|')
                isfname      = att[0].strip()
                clr[isfname] = att[1].strip()
                lstyle=True
        if not lstyle:
            print runid+' not found in isf.sty'
            raise Exception
    except Exception as e:
        print e
        sys.exit(42)

    return clr


# check argument
def load_argument():
    # deals with argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-runid", metavar='runid list' , help="used to look information in runid.db"                  , type=str, nargs='+' , required=True )
    parser.add_argument("-f"    , metavar='file list'  , help="file list to plot (default is runid_var.nc)"           , type=str, nargs='+' , required=False)
    parser.add_argument("-var"  , metavar='var list'   , help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs='+' , required=True )
    parser.add_argument("-varf" , metavar='file list'  , help="file list to plot (1 file pattern per variable)"       , type=str, nargs='+' , required=False)
    parser.add_argument("-title", metavar='title'      , help="subplot title (associated with var)"                   , type=str, nargs='+' , required=False)
    parser.add_argument("-dir"  , metavar='directory of input file' , help="directory of input file"                  , type=str, nargs=1   , required=False, default=['./'])
    parser.add_argument("-sf"  , metavar='scale factor', help="scale factor"                                          , type=float, nargs=1 , required=False, default=[1])
    parser.add_argument("-minmax", metavar='min / max' , help="min/max"                                               , type=float, nargs=2 , required=False, default=[1])
    parser.add_argument("-o"    , metavar='figure_name', help="output figure name without extension"                  , type=str, nargs=1   , required=False, default=['output'])
    # flag argument
    parser.add_argument("-obs"  , metavar='obs mean and std file', help="obs mean and std file"          , type=str, nargs='+', required=False)
    parser.add_argument("-mean" , help="will plot model mean base on input netcdf file"                               , required=False, action="store_true")
    parser.add_argument("-noshow" , help="do not display the figure (only save it)"                                   , required=False, action="store_true")
    return parser.parse_args()

def output_argument_lst(cfile, arglst):
    fid = open(cfile, "w")
    fid.write(' python2.7 '+' '.join(arglst))
    fid.close()

# ============================ plotting tools ==================================
def get_corner(ax):
    x0 = ax.get_position().x1
    x1 = x0+0.1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1
    return x0, x1, y0, y1

def add_legend(lg, ax, lvar, lrun, ncol=7, lvis=True):
    x0, x1, y0, y1 = get_corner(ax)
    lline, llabel = lg.get_legend_handles_labels()
    # isf name
    lax = plt.axes([0.00, 0.07, 1, 0.08])
    leg=plt.legend(lline[0::len(lrun)], lvar, loc='upper left', ncol = ncol, fontsize=16, frameon=False)
    for item in leg.legendHandles:
        item.set_visible(lvis)
    lax.set_axis_off() 
    # run name
    lax = plt.axes([0.00, 0.0, 1, 0.06])
    leg=plt.legend(lline[0:len(lrun)], lrun, loc='upper left', ncol = ncol, fontsize=16, frameon=False)
    for item in leg.legendHandles:
        item.set_visible(lvis)
    lax.set_axis_off() 

def add_text(lg, ax, cisf, crun, ncol=3, lvis=True):
    x0, x1, y0, y1 = get_corner(ax)
    lline, llabel = lg.get_legend_handles_labels()
    lax = plt.axes([0.07, 0.0, 1, 0.15])
    leg=plt.legend(lline[0:end:nrun], cisf, loc='upper left', ncol = ncol, fontsize=16, frameon=False)
    for item in leg.legendHandles:
        item.set_visible(lvis)
    lax.set_axis_off() 

    lax = plt.axes([0.0, 0.0, 1, 0.06])
    leg=plt.legend(lline[0:nrun], crun, loc='upper left', ncol = ncol, fontsize=16, frameon=False)
    for item in leg.legendHandles:
        item.set_visible(lvis)
    lax.set_axis_off() 
# ========================== stat plot ============================
def tidyup_ax(ax, xmin, xmax, ymin, ymax):
    ax.set_ylim([ymin, ymax])
    ax.set_xlim([xmin, xmax])
    ax.grid(visible=True)
    ax.set_yticklabels([])
    ax.set_xticks([])

#def add_modstat(ax, run_lst):
#    for irun in range(len(run_lst)):
#        cpl = plt.errorbar(irun+1, run_lst[irun].mean, yerr=run_lst[irun].std, fmt='o', markeredgecolor=run_lst[irun].color, markersize=8, color=run_lst[irun].color, linewidth=2)

def add_obsstat(ax, cobs, cisf, x_isf):
    mean = cobs.dta[cisf][0] ; std = cobs.dta[cisf][1]
    cpl = plt.errorbar(x_isf, mean, yerr=std, fmt=cobs.mkr, markeredgecolor=cobs.clr[cisf], markersize=8, color=cobs.clr[cisf], linewidth=2)

# ============================ file parser =====================================
def parse_dbfile(runid, cisf):
    try:
        lstyle=False
        with open('run.sty') as fid:
            for cline in fid:
                att=cline.split('|')
                if att[0].strip() == runid:
                    cpltrunid = att[0].strip()
                    cpltname  = att[1].strip()
                    cpltline  = att[2].strip()
                    lstyle=True
        if not lstyle:
            print runid+' not found in run.sty'
            raise Exception
       
        lstyle=False
        with open('isf.sty') as fid:
            for cline in fid:
                att=cline.split('|')
                if att[0].strip() == cisf:
                    cpltcolor  = att[1].strip()
                    lstyle=True
        if not lstyle:
            print runid+' not found in isf.sty'
            raise Exception


    except Exception as e:
        print e
        sys.exit(42)

    # return value
    return cpltrunid, cpltname, cpltline, cpltcolor
# ============================ file parser end =====================================

def main():

# load argument
    args = load_argument()

# output argument list
    output_argument_lst(args.o[0]+'.txt', sys.argv)

# initialisation
    nrun = len(args.runid)
    nvar = len(args.var)
    crun = [None]*nrun
    nobs = len(args.obs)
    obs_lst  = [None]*nobs
    ax       = [None]
    rmin = args.minmax[0] ; rmax = args.minmax[1]
    disf=load_isfclr()

    plt.figure(figsize=np.array([290, 210]) / 25.4)
 
# need to deal with multivar
    mintime=dt.date.max
    maxtime=dt.date.min
    ymin=-sys.float_info.max
    ymax=sys.float_info.max

    ax[0] = plt.subplot(1, 1, 1)
    if args.obs:
        for iobs, cobs in enumerate(args.obs):
        # load obs
            obs_lst[iobs] = obs(cobs, disf)   # name, mkr, dta
 
    for ivar, cvar in enumerate(args.var):
        for irun, runid in enumerate(args.runid):

            rundta = run(runid, cvar, args.sf[0])
            crun[irun] = rundta.name

            # load data
            if args.f:
                if len(args.f) == 1 :
                    fglob = args.f[0]
                else :
                    fglob = args.f[irun]
                cfile = glob.glob(args.dir[0]+'/'+runid+'/'+fglob)
                if len(cfile)==0:
                    print 'no file found with this pattern '+args.dir[0]+'/'+runid+'/'+fglob
                    sys.exit(42)
            elif args.varf:
                if len(args.varf) == 1 :
                    fglob = args.varf[0]
                else :
                    fglob = args.varf[ivar]
                cfile = glob.glob(args.dir[0]+'/'+runid+'/'+fglob)
                if len(cfile)==0:
                    print 'no file found with this pattern '+args.dir[0]+'/'+runid+'/'+fglob
                    sys.exit(42)
            else:
                cfile = glob.glob(args.dir[0]+'/'+runid+'_'+cvar+'.nc')
                if len(cfile)==0:
                    print 'no file found with this pattern '+args.dir[0]+'/'+runid+'_'+cvar+'.nc'
                    sys.exit(42)

            rundta.load_time_series(cfile, 'isfmelt_'+cvar)
            ts = rundta.ts
            lg = ts.plot(ax=ax[0], legend=False, style=rundta.line, color=disf[cvar], markeredgecolor=disf[cvar], markersize=8, x_compat=True, linewidth=2, rot=0)
            #
            # limit of time axis
            mintime=min([mintime,ts.index[0].to_pydatetime().date()])
            maxtime=max([maxtime,ts.index[-1].to_pydatetime().date()])

        # set title
        if (args.title):
            ax[0].set_title(args.title[0],fontsize=20)

        # set x axis
        nlabel=5
        ndays=(maxtime-mintime).days
        nyear=ndays/365
        if nyear < 10:
            nyt=1
        elif 10<=nyear<50:
            nyt=5
        elif 50<=nyear<100:
            nyt=10
        else:
            nyt=100

        ax[0].xaxis.set_major_locator(mdates.YearLocator(nyt,month=1,day=1))
        ax[0].tick_params(axis='both', labelsize=16)
        if (ivar != nvar-1):
            ax[0].set_xticklabels([])
        else:
            ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

        for lt in ax[0].get_xticklabels():
            lt.set_ha('center')

    ax[0].set_ylim([rmin,rmax])
    ax[0].grid()
 
    # tidy up space
    plt.subplots_adjust(left=0.1, right=0.8, bottom=0.2, top=0.92, wspace=0.15, hspace=0.15)

    # add legend
    add_legend(lg,ax[0],args.var, crun)
    print lg

    if args.mean or args.obs:
        xmin = -1./nvar ; xmax = 1.+1./nvar
        for ivar, cvar in enumerate(args.var):
            x0, x1, y0, y1=get_corner(ax[0])    
            cax = plt.axes([x0+0.01, y0, x1-x0, y1-y0])
            # set min/max/grid ...
            tidyup_ax(cax, xmin, xmax, rmin, rmax)

            # plot obs mean
            if args.obs:
                for iobs, cobs in enumerate(args.obs):
                    add_obsstat(cax, obs_lst[iobs], cvar, 1.*ivar/nvar)

    plt.savefig(args.o[0]+'.png', format='png', dpi=150)

    if args.noshow: 
       pass
    else:
       plt.show()

if __name__=="__main__":
    main()
