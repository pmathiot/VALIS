import numpy as np
import math
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
import cartopy
import cartopy.crs as ccrs

# def class runid
#GETZ_CTD_WOA2018.nc
class obs(object):
    def __init__(self,isfc,mshc):
        self.name='WOA2018'
        self.cfile=isfc.name+'_CTD_WOA2018.nc'
        self.cast, self.lon, self.lat, self.ncast=self.load_cast(isfc,mshc)
        self.obsT,self.z=self.load_prof('Temperature','z')
        self.obsm=self.load_mltobs(isfc)

    def load_mltobs(self,isfc):
        with open('Rignot_2013.txt_test', 'r') as obsfile:
            for line in obsfile:
                isfname    = line.split()[0]
                if isfname == isfc.name:
                    isfmlt_val = float(line.split()[1])
                    isfmlt_err = float(line.split()[2])
                    break
        return (isfmlt_val, isfmlt_err)
        
    def load_cast(self,isfc,mshc):
        valid_cast=[]
        valid_lon=[]
        valid_lat=[]
        ncid = nc.Dataset(self.cfile)
        lon = ncid.variables['lon'][:]
        lat = ncid.variables['lat'][:]
        flgp = ncid.variables['Temperature_WODprofileflag'][:]
        ncid.close()
        ncast = len(lat)
        for icast in range(0,ncast):
            if ( flgp[icast] == 0 ):
                d,jmod,imod = get_shortest_in(lon[icast], lat[icast], mshc.lon, mshc.lat, isfc.msk+isfc.msk_cnt)
                if (d < 30):
                    valid_lon.append(lon[icast])
                    valid_lat.append(lat[icast])
                    valid_cast.append(icast)
                    print lon[icast], lat[icast], icast, d
        return valid_cast,valid_lon,valid_lat,len(valid_cast)

    def load_prof(self,cvarT,cvarz):
        ncid = nc.Dataset(self.cfile)
        cast_ndta = ncid.variables['Temperature_row_size'][:]
        Temp_dta  = ncid.variables['Temperature'][:].data
        z_dta  = ncid.variables['z'][:]
        flg = ncid.variables['Temperature_WODflag'][:]
        profz=np.empty(shape=(136,self.ncast)); profz.fill(np.nan)
        profT=np.empty(shape=(136,self.ncast)); profT.fill(np.nan)
        flgc=np.empty(shape=(len(flg))); flgc.fill(np.nan)
        flgc[(flg==0) & (np.abs(Temp_dta) < 100)]=1
        for iprof, icast in enumerate(self.cast):
            ikmin=0 ; ikmax=cast_ndta[icast]
            if icast == 0:
                iimin=0
            else:
                iimin=np.sum(cast_ndta[:icast])
            iimax=iimin+ikmax
            profz[ikmin:ikmax,iprof]=z_dta[iimin:iimax]*flgc[iimin:iimax]
            profT[ikmin:ikmax,iprof]=Temp_dta[iimin:iimax]*flgc[iimin:iimax]
        return profT,profz

# compute distance
def get_shortest_in(xlon, ylat, lon2d, lat2d, msk2d):
    """needle is a single (lat,long) tuple.
        haystack is a numpy array to find the point in
        that has the shortest distance to needle
    """
    earth_radius_miles = 3956.0
    dlat = np.radians(lat2d[:,:]) - np.radians(ylat)
    dlon = np.radians(lon2d[:,:]) - np.radians(xlon)
    a = np.square(np.sin(dlat/2.0)) + np.cos(np.radians(ylat)) * np.cos(np.radians(lat2d[:,:])) * np.square(np.sin(dlon/2.0))
    [jmin, imin] = np.nonzero(a==np.amin(a[msk2d==1]))   # the minimun location can be retreive at this stage
    great_circle_distance = 2 * np.arcsin(np.sqrt(a[jmin, imin]))
    d = earth_radius_miles * great_circle_distance
    return d, jmin, imin

class run(object):
    def __init__(self, runid, cdir, cfile, cfilets, cvmlt, cvtem, isfc, mshc):
        # parse dbfile
        self.dnc   = {}
        self.isfd  = {}
        self.runid = runid
        self.name, self.clr  = self.load_att(runid)

        self.load_dnc(cdir, cfile, cvmlt, '1', 'mltd')
        self.load_dnc(cdir, cfile, cvtem, '1', 'temp')
        self.load_dnc(cdir, cfilets, cvmlt, '1', 'mltts')
        self.draft, self.area, self.melt, self.T, self.meltts, self.melttot = self.load_data(isfc,mshc)

    def load_dnc(self, cdir, cfile, cvar, cnf, ckey):
        cfilel = sorted(glob.glob(cdir+'/'+self.runid+'/'+cfile))
        if len(cfilel)==0:
            print 'no file found with this pattern '+cdir+'/'+self.runid+'/'+cfile
            sys.exit(42)

        self.dnc['cf'+ckey]=cfilel
        self.dnc['cv'+ckey]=cvar

    def load_data(self, isfc, mshc):
        area = mshc.area[isfc.msk==1]
        scal = 86400. * 365. / (1000. * 1e9) * -1.

        # mlt distri
        drft = mshc.zisf[isfc.msk==1]
        mltr = get_2d_data(self.dnc['cfmltd'][0],self.dnc['cvmltd'])
        mltr = np.ma.masked_where(mltr == 0.0, mltr).filled(np.nan) 
        mltd = mltr[isfc.msk==1] * area * scal

        # mlt total
        mltt = np.sum(mltr[(isfc.msk==1)] * mshc.area[(isfc.msk==1)] * scal)
       
        # temp profile
        nk = get_dim(self.dnc['cftemp'][0],'z')
        proT=np.zeros(shape=(nk,2))
        for jk in range(0,nk):
            temp2d = get_2d_data(self.dnc['cftemp'][0],self.dnc['cvtemp'],klvl=jk) 
            msk2d  = get_2d_data('mesh_mask.nc','tmask',klvl=jk) 
            proT[jk,0] = mshc.z[jk] 
            proT[jk,1] = np.ma.average(temp2d,weights=mshc.area*isfc.msk_cnt*msk2d)

        # melt ts
        mltts=[] ; timts=[]
        for cfile in self.dnc['cfmltts']:
            timts.append(get_time_data(cfile,'time_centered')[:])
            mltr = get_2d_data(cfile,self.dnc['cvmltts'])
            mltts.append(np.sum(mltr[(isfc.msk==1)] * mshc.area[(isfc.msk==1)] * scal ))

        return drft,area/(1000*1000),mltd,proT, [timts, mltts], mltt

    def load_att(self,runid):
        try:
            lstyle=False
            with open('run.sty') as fid:
                for cline in fid:
                    att=cline.split('|')
                    if att[0].strip() == runid:
                        cpltname  = att[1].strip()
                        cpltclr   = att[3].strip()
                        lstyle=True
            if not lstyle:
                print runid+' not found in run.sty'
                raise Exception
           
        except Exception as e:
            print 'load_name error on ',runid
            print e
            sys.exit(42)
    
        # return value
        return cpltname, cpltclr

    def __str__(self):
        return 'runid = {}, name = {}'.format(self.runid, self.name)

#================================================================================
class isf(object):
    def __init__(self, cisf):
        self.name = cisf
        self.clr, self.iid, self.mapdef = self.load_isf(cisf)
        self.msk     = self.load_msk(cisf,'mask_isf.nc','mask_isf')
        self.msk_cnt = self.load_msk(cisf,'mask_isf.nc','mask_isf_front')

    def load_msk(self, cisf, cfile, cvar):
        lmsk=get_2d_data(cfile,cvar)
        lmsk[lmsk != -self.iid] = 0
        lmsk[lmsk == -self.iid] = 1
        return lmsk

    def load_isf(self, cisf):
        try:
            lstyle=False
            with open('isf.sty') as fid:
                for cline in fid:
                    att=cline.split('|')
                    isfname      = att[0].strip()
                    if isfname == cisf :
                        clr = att[1].strip()
                        iid = int(att[2].strip())
                        mapdef = [np.asarray(eval(att[3].strip())),ccrs.Stereographic(central_latitude=-90.0, central_longitude=0)]
                    lstyle=True
            if not lstyle:
                print runid+' not found in isf.sty'
                raise Exception
    
            return clr,iid,mapdef

        except Exception as e:
            print e
            sys.exit(42)

#================================================================================
class msh(object):
    def __init__(self):
        ncid   = nc.Dataset('mesh_mask.nc')
        self.z    = ncid.variables['gdept_1d'][:].squeeze()
        ncid.close()

        self.lon  = get_2d_data('coordinates.nc','glamt')
        self.lat  = get_2d_data('coordinates.nc','gphit')

        self.zisf = get_2d_data('mesh_mask.nc','isfdraft')
        self.area = get_2d_data('mesh_mask.nc','e1t') * get_2d_data('mesh_mask.nc','e2t')

#================================================================================
class plot(object):
    def __init__(self, ctitle, cfout,isfc):
        self.ax={}
        self.title=ctitle
        self.fout=cfout
        self.fig = plt.figure(figsize=np.array([290, 210]) / 25.4)
        self.ax['aread'] = self.fig.add_subplot(1, 2, 2)
        self.ax['meltd'] = self.ax['aread'].twiny()
        self.ax['profT'] = self.fig.add_subplot(1, 2, 1)
        plt.subplots_adjust(left=0.1, right=0.82, bottom=0.52, top=0.88, wspace=0.15, hspace=0.15)
        self.ax['meltt'] = plt.axes([0.92, 0.52, 0.06 , 0.36])
        self.ax['meltts']= plt.axes([0.1 , 0.10, 0.56 , 0.34])
        self.ax['mapisf']= plt.axes([0.68, 0.10, 0.30 , 0.34],projection=isfc.mapdef[1])

    def plot_dist(self,runlst,isfc,obsc,mshc):

        self.plot_aread(runlst)
        self.plot_meltd(runlst)
        self.plot_meltt(runlst,obsc)
        self.plot_meltts(runlst,obsc)
        self.plot_obst(obsc)
        self.plot_proft(runlst)
        self.plot_mapisf(isfc,obsc,mshc)

        self.add_legend()

        self.add_title()

        self.set_axes('aread',False,'km2','top', ylim=[0, 1000])
        self.set_axes('meltd',True,'Gt/y','bottom', ylim=[0, 1000])
        self.set_axes('profT',True,'C','bottom',cylabel='Depth [m]', ylim=[0, 1000])
        self.set_axes('meltt',False,'Gt/y','bottom',cylabel='Total melt [Gt/y]', lyrevert=False)
        self.set_axes('meltts',True,'Gt/y','bottom',cylabel='Total melt [Gt/y]', lyrevert=False)
        self.set_axes('mapisf',False,'','bottom', lyrevert=False)

        self.save()

    def plot_mapisf(self,isfc,obsc,mshc):
        isf_features   = cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines', '50m',facecolor='none',edgecolor='k')
        coast_features = cartopy.feature.NaturalEarthFeature('physical', 'coastline'                         , '50m',facecolor='0.75',edgecolor='k')
        bathy1000_features = cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_J_1000'          , '10m',facecolor='none',edgecolor='k')
        bathy2000_features = cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_I_2000'          , '10m',facecolor='none',edgecolor='k')
        bathy3000_features = cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_H_3000'          , '10m',facecolor='none',edgecolor='k')

        cmap=plt.get_cmap('Oranges',2)
        cmap.set_under('w', 1.0)
        coords = isfc.mapdef[1].transform_points( ccrs.PlateCarree(), isfc.mapdef[0][0:2], isfc.mapdef[0][2:4])
        self.ax['mapisf'].set_extent([coords[0, 0], coords[1, 0], coords[0, 1], coords[1, 1]], isfc.mapdef[1])
        self.ax['mapisf'].gridlines()
        for feature in [coast_features, isf_features]:#, bathy1000_features, bathy2000_features, bathy3000_features]:
            self.ax['mapisf'].add_feature(feature)
        self.ax['mapisf'].pcolormesh(mshc.lon[0:-4,:],mshc.lat[0:-4,:],isfc.msk[0:-4,:],cmap=cmap,vmin=0.5, vmax=1.5,transform=ccrs.PlateCarree(),rasterized=True, label=isfc.name)
        self.ax['mapisf'].plot(obsc.lon,obsc.lat,'o', markersize=6, color='darkturquoise', markeredgecolor='darkturquoise', transform=ccrs.PlateCarree(), label='WOA2018 obs')

    def plot_proft(self,runlst):
        for irun, runid in enumerate(runlst):
            self.ax['profT'].plot(runid.T[:,1],runid.T[:,0],color=runid.clr,linewidth=2,label=runid.name)

    def plot_obst(self,obsc):
        for icast in range(0,obsc.ncast):
            self.ax['profT'].plot(obsc.obsT[:,icast],obsc.z[:,icast],color='gainsboro',linewidth=2)
        self.ax['profT'].plot(np.nanmean(obsc.obsT,axis=1),np.nanmean(obsc.z,axis=1),color='silver',linewidth=5,label=obsc.name)

    def plot_aread(self,runlst):
        # plot area distri only for one simulation
        runid = runlst[0]
        self.ax['aread'].hist(runid.draft, weights=runid.area, bins=20, range=(0,1000),   \
        orientation='horizontal',color='lightgray', edgecolor='lightgray',linewidth=2, label='isf area')

    def plot_meltd(self,runlst):
        for irun, runid in enumerate(runlst):
            self.ax['meltd'].hist(runid.draft, weights=runid.melt, bins=20, range=(0,1000),   \
            orientation='horizontal',facecolor="None", edgecolor=runid.clr,linewidth=2)

    def plot_meltt(self,runlst,obsc):
        self.ax['meltt'].errorbar(0.,obsc.obsm[0], yerr=obsc.obsm[1], fmt='o', markersize=8, color='k', linewidth=3)
        for irun, runid in enumerate(runlst):
            self.ax['meltt'].plot(runid.melttot, marker='o', markeredgecolor=runid.clr, color=runid.clr, markersize=8)

    def plot_meltts(self,runlst,obsc):
        timemax=np.max([runlst[irun].meltts[0][-2] for irun in range(0,len(runlst))])
        for irun, runid in enumerate(runlst):
            self.ax['meltts'].plot(runid.meltts[0][:], runid.meltts[1][:], '-o', markeredgecolor=runid.clr, color=runid.clr, markersize=8)
        self.ax['meltts'].errorbar(timemax, obsc.obsm[0], yerr=obsc.obsm[1], fmt='o', markersize=8, color='k', linewidth=3, label='Rignot et al. 2013')

    def add_title(self):
        self.ax['profT'].set_title('Temperature (isf front) [C]',fontsize=18)
        self.ax['meltd'].set_title('Ice shelf melt [Gt/y]',fontsize=18)
        self.fig.suptitle(self.title,fontsize=24)

    # set axes bound and grid
    def set_axes(self,cname,lxtick,cxlabel,cxpos, cylabel=None, lyrevert=True, ylim=None):
        if (ylim):
           self.ax[cname].set_ylim(ylim)
        if (lyrevert):
           self.ax[cname].invert_yaxis()
        self.ax[cname].tick_params(labelsize=16)
        self.ax[cname].xaxis.set_ticks_position(cxpos)
        #self.ax[cname].xaxis.set_label_position(cxpos)
        self.ax[cname].grid()
        self.ax[cname].grid(visible=True)
        #self.ax[cname].set_xlabel(cxlabel, fontsize=16)

        if not lxtick:
            self.ax[cname].set_xticks([])

        if cylabel:
            self.ax[cname].set_ylabel(cylabel, fontsize=16)
        else:
            self.ax[cname].set_yticklabels('')

    # tidy up space
        #plt.subplots_adjust(left=0.1, right=0.95, bottom=0.12, top=0.88, wspace=0.15, hspace=0.15)

    # save plot
    def save(self):
        self.fig.savefig(self.fout+'.png', format='png', dpi=150)

    def add_legend(self, ncol=5, lvis=True):
        # isf area
        lline, llabel = self.ax['aread'].get_legend_handles_labels()
        leg = self.ax['aread'].legend(lline, llabel, loc='upper right', fontsize=16, frameon=False)
        for item in leg.legendHandles:
            item.set_visible(lvis)

        # profile obs
        lline, llabel = self.ax['profT'].get_legend_handles_labels()
        leg=  self.ax['profT'].legend([lline[0]], [llabel[0]], loc='upper right', fontsize=16, frameon=False)

        # isf melt
        lline, llabel = self.ax['meltts'].get_legend_handles_labels()
        leg=  self.ax['meltts'].legend([lline[-1]], [llabel[-1]], numpoints=1, loc='upper left', fontsize=16, frameon=False)

        # isf melt
        lline, llabel = self.ax['mapisf'].get_legend_handles_labels()
        leg=  self.ax['mapisf'].legend(lline, llabel, numpoints=1, loc='upper left', fontsize=16, frameon=False)

        lline, llabel = self.ax['profT'].get_legend_handles_labels()
        lax = plt.axes([0.00, 0.00, 1, 0.08])
        leg=  plt.legend(lline[1::], llabel[1::], loc='upper left', ncol = ncol, fontsize=16, frameon=False)
        for item in leg.legendHandles:
            item.set_visible(lvis)
        lax.set_axis_off()

# ================================= NETCDF ===============================
def get_name(regex,varlst):
    revar = re.compile(r'\b%s\b'%regex,re.I)
    cvar  = filter(revar.match, varlst)
    if (len(cvar) > 1):
        print regex+' name list is longer than 1 or 0; error'
        print cvar
        print cvar[0]+' is selected'
    if (len(cvar) == 0):
        print 'no match between '+regex+' and :'
        print varlst
        sys.exit(42)
    return cvar[0]

def get_variable_shape(ncid,ncvar):
    redimt=re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    redimz=re.compile(r"\b(z|dep|nav_lev|depth|deptht)\b", re.I)
    redimy=re.compile(r"\b(y|y_grid_.+|latitude|lat|nj)\b", re.I)
    redimx=re.compile(r"\b(x|x_grid_.+|lon|longitude|long|ni)\b", re.I)
    dimlst = ncvar.dimensions
    if (len(ncvar.shape)==2) and redimx.match(dimlst[1]) and redimy.match(dimlst[0]):
        cshape='XY'
    elif (len(ncvar.shape)==3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape='XYT'
    elif (len(ncvar.shape)==3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimz.match(dimlst[0]):
        cshape='XYZ'
    elif (len(ncvar.shape)==4) and redimx.match(dimlst[3]) and redimy.match(dimlst[2]) and redimz.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape='XYZT'
    else:
        print 'cshape undefined, error'
        print dimlst
        sys.exit(42)
    return cshape

def get_dim(cfile,cdir):
    ncid   = nc.Dataset(cfile)
    if cdir=='x' :
        redim=re.compile(r"\b(x|x_grid_.+|lon|longitude|long)\b", re.I)
    elif cdir=='y' :
        redim=re.compile(r"\b(y|y_grid_.+|latitude|lat)\b", re.I)
    elif cdir=='z' :
        redim=re.compile(r"\b(z|dep|nav_lev|depth|deptht)\b", re.I)
    elif cdir=='t' :
        redim=re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    else:
        print 'dimension direction unknown, need to be x, y, z or k'
        sys.exit(42)

    cdim=filter(redim.match, ncid.dimensions.keys());
    if (len(cdim) > 1):
        print regex+' name list is longer than 1; error'
        print cdim
        sys.exit(42)
    elif (len(cdim) == 0):
        print cdir+' dim in '+cfile+' is 0.'
        ndim=0
    else:
        cdim=cdim[0]
        ndim=len(ncid.dimensions[cdim])
    ncid.close()
    return ndim

def get_dims(cfile):
    nx = get_dim(cfile,'x')
    ny = get_dim(cfile,'y')
    nz = get_dim(cfile,'z')
    nt = get_dim(cfile,'t')
    return nx,ny,nz,nt

def get_time_data(cfile,cvar):
    ncid    = nc.Dataset(cfile)

    ncvtime = ncid.variables[cvar]
    if 'units' in ncvtime.ncattrs():
        cunits = ncvtime.units
    else:
        cunits = "seconds since 1900-01-01 00:00:00"

    # define calendar
    if 'calendar' in ncvtime.ncattrs():
        ccalendar = ncvtime.calendar
    else:
        ccalendar = "noleap"
    time = nc.num2date(ncid.variables[cvar][:].squeeze(), cunits, ccalendar)

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

    return timeidx[:]

# get_2d_data
def get_2d_data(cfile,cvar,ktime=0,klvl=0):
    print ' reading '+cvar+' in '+cfile+' ...'

    nx,ny,nz,nt=get_dims(cfile)

    ncid  = nc.Dataset(cfile)
    clvar = get_name(cvar,ncid.variables.keys())
    var   = ncid.variables[clvar]
    shape = get_variable_shape(ncid,var)

    if shape=='XY' :
        print ' 2d variable XY'
        if (klvl > 0) :
            print 'error klvl larger than 0 (klvl = '+str(klvl)+')'
            sys.exit(42)
        if (ktime > 0) :
            print 'error ktime larger than 0 (ktime = '+str(ktime)+')'
            sys.exit(42)
        dat2d=var[:,:]
    elif shape=='XYT' :
        print ' 3d variable XYT'
        if (klvl > 0) :
            print 'error klvl larger than 0 (klvl = '+str(klvl)+')'
            sys.exit(42)
        if (ktime > 0) :
            print 'error ktime larger than 0 (ktime = '+str(ktime)+')'
            sys.exit(42)
        dat2d=var[ktime,:,:]
    elif shape=='XYZ' :
        print ' 3d variable XYZ'
        if (ktime > 0) :
            print 'error ktime larger than 0 (ktime = '+str(ktime)+')'
            sys.exit(42)
        dat2d=var[klvl,:,:]
    elif len(var.shape)==4 :
        print ' 4d variable XYZT'
        dat2d=var[ktime,klvl,:,:]
    else:
        print cvar+' contains '+str(len(var.shape))+' dimensions'
        print 'dimension names are '+var.dimensions
        print ' shape unknown, exit '
        sys.exit(42)
    ncid.close()
    return dat2d

#================================================================================
# check argument
def load_argument():
    # deals with argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-runid", metavar='runid list' , help="used to look information in runid.db"                  , type=str, nargs="+" , required=True )
    parser.add_argument("-fmlt" , metavar='file list'  , help="file list to plot (default is runid_var.nc)"           , type=str, nargs="+" , required=False)
    parser.add_argument("-fmltts", metavar='file list ts'  , help="file list to plot (default is runid_var.nc)"       , type=str, nargs="+" , required=False)
    parser.add_argument("-vmlt" , metavar='isf mlt var', help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs="+" , required=True )
    parser.add_argument("-vtem" , metavar='oce tem var', help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs="+" , required=True )
    parser.add_argument("-isf"  , metavar='isf name'   , help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs=1   , required=True )
    parser.add_argument("-dir"  , metavar='directory of input file' , help="directory of input file"                  , type=str, nargs=1   , required=False, default=['./'])
    parser.add_argument("-title", metavar='figure title', help="figure title"                                         , type=str, nargs=1   , required=True)
    parser.add_argument("-o"    , metavar='figure_name', help="output figure name without extension"                  , type=str, nargs=1   , required=False, default=['output'])
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

def fix_arglist(nrun,argin,ctxt):
    argout=[None]*nrun
    if (len(argin) == 1):
       argout = argin * nrun
       print 'all argument ('+ctxt+') set to the same value'
    elif (len(argin) == nrun):
       argout = argin
    else:
       print 'number of arguments is not compatible with the number of runid ('+ctxt+')'
       sys.exit(42)
    return argout


def main():

# load argument
    args = load_argument()

# output argument list
    output_argument_lst(args.o[0]+'.txt', sys.argv)

# load msh data
    mshc = msh()

# load isf data
    isfc = isf(args.isf[0])

# load WOA profile
    profc = obs(isfc,mshc)

# initialisation
    plt_isf = plot(args.title[0], args.o[0],isfc)

# load data for each run
    nrun=len(args.runid)
    run_lst=[None]*nrun

    fmlt  =fix_arglist(nrun,args.fmlt  ,'args.fmlt')
    fmltts=fix_arglist(nrun,args.fmltts,'args.fmltts')
    vmlt  =fix_arglist(nrun,args.vmlt  ,'args.vmlt')
    vtem  =fix_arglist(nrun,args.vtem  ,'args.vtem')
    cdir  =fix_arglist(nrun,args.dir   ,'args.dir')

    for irun, runid in enumerate(args.runid):
        run_lst[irun] = run(runid, cdir[irun], fmlt[irun], fmltts[irun], vmlt[irun], vtem[irun], isfc, mshc)

    plt_isf.plot_dist(run_lst,isfc,profc,mshc)

    if args.noshow:
       pass
    else:
       plt.show()

if __name__=="__main__":
    main()
