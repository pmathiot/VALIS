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
import gsw

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

# def class runid
#===========================================================================
class obs(object):
    def __init__(self,isfc,mshc):
        self.name='WOA2018'
        self.cfile=isfc.name+'_CTD_WOA2018.nc'

        self.cast, self.lon, self.lat, self.ncast = self.load_cast('Temperature',isfc,mshc)
        self.obsT,self.Tz=self.load_prof('Temperature')
        self.Tncast=self.ncast
        self.Tcast =self.cast

        self.cast, self.lon, self.lat, self.ncast = self.load_cast('Salinity',isfc,mshc)
        self.obsS,self.Sz=self.load_prof('Salinity')
        self.Sncast=self.ncast
        self.Scast =self.cast

        self.obsm=self.load_mltobs(isfc)

    def load_mltobs(self,isfc):
        with open('MLT/Rignot_2013.txt', 'r') as obsfile:
            for line in obsfile:
                isfname    = line.split()[0]
                if isfname == isfc.name:
                    isfmlt_val = float(line.split()[1])
                    isfmlt_err = float(line.split()[2])
                    break
        return (isfmlt_val, isfmlt_err)

    def load_prof(self,cvar):
        #
        # load data
        ncid      = nc.Dataset('OBS/'+self.cfile)

        cast_vndta= ncid.variables[cvar+'_row_size'][:]
        flg       = ncid.variables[cvar+'_WODflag'][:]
        var_dta   = ncid.variables[cvar][:].data

        cast_zndta= ncid.variables['z_row_size'][:]
        z_dta     = ncid.variables['z'][:]
        ncid.close()
        #
        # fill mask array
        flgc = np.where(flg==0,1.,np.nan)
        #
        # fill profile
        profz=np.empty(shape=(136,self.ncast)); profz.fill(np.nan)
        profv=np.empty(shape=(136,self.ncast)); profv.fill(np.nan)
        for ic, cast_id in enumerate(self.cast):
            ikvmax=cast_vndta[cast_id]
            ikzmax=cast_zndta[cast_id]
            iivmin=np.sum(cast_vndta[:cast_id], dtype=np.int32) # if icast = 0 => iimin=0 by np.sum properties
            iizmin=np.sum(cast_zndta[:cast_id], dtype=np.int32) # if icast = 0 => iimin=0 by np.sum properties
            slice_vdta=slice(iivmin,iivmin+ikvmax)
            slice_zdta=slice(iizmin,iizmin+ikzmax)
            profv[:ikvmax,ic]=var_dta[slice_vdta]*flgc[slice_vdta]
            profz[:ikzmax,ic]=z_dta  [slice_zdta]#*flgc[slice_zdta]
        return profv,profz
       
    def load_cast(self,cvar,isfc,mshc):
        #
        # load data set
        ncid      = nc.Dataset('OBS/'+self.cfile)
        lon       = ncid.variables['lon'][:]
        lat       = ncid.variables['lat'][:]
        flgp      = ncid.variables[cvar+'_WODprofileflag'][:]
        cast_size = ncid.variables[cvar+'_row_size'][:]
        ncid.close()
        #
        # only work on valid profile
        valid_flag_idx=np.where( (flgp==0) & (cast_size>0))
        #
        # compute distance (very slow)
        d=np.ones(shape=len(lat))*1.e20
        msk_idx=np.where(isfc.msk_cnt*mshc.msk==1)
        mshc_lonr=np.radians(mshc.lon[msk_idx])
        mshc_latr=np.radians(mshc.lat[msk_idx])
        lonr = np.radians(lon)
        latr = np.radians(lat)
        for icast in valid_flag_idx[0]:
            d[icast]=get_shortest_in(lonr[icast], latr[icast], mshc_lonr, mshc_latr)
        #
        # extract location and cast number
        # only data close to the isf
        valid_lon =[]
        valid_lat =[]
        valid_dist_idx=np.where(d<30)
        valid_lon.append(lon[valid_dist_idx])
        valid_lat.append(lat[valid_dist_idx])
        valid_cast=valid_dist_idx[0]

        return valid_cast,valid_lon,valid_lat,len(valid_cast)

# compute distance
def get_shortest_in(xlon, ylat, lon1d, lat1d):
    """needle is a single (lat,long) tuple.
        haystack is a numpy array to find the point in
        that has the shortest distance to needle
    """
    dlat = lat1d - ylat
    dlon = lon1d - xlon
    a = np.square(np.sin(dlat/2.0)) + np.cos(ylat) * np.cos(lat1d) * np.square(np.sin(dlon/2.0))
    d = 6371.0 * 2.0 * np.arcsin(np.sqrt(a[a==np.amin(a)]))
    return d

#===========================================================================
class run(object):
    def __init__(self, runid, cdir, cfiletts, cfilests, cfileqts, cvtem, cvsal, cvmlt, isfc, mshc):
        # parse dbfile
        self.cdir  = cdir
        self.dnc   = {}
        self.isfd  = {}
        self.runid = runid
        self.name, self.clr  = self.load_att(runid)

        self.load_dnc(cdir, cfiletts, cvtem, '1', 'tempts')
        self.load_dnc(cdir, cfilests, cvsal, '1', 'salts' )
        self.load_dnc(cdir, cfileqts, cvmlt, '1', 'mltts' )
        self.draft, self.area, self.melt, self.meltts, self.melttot, self.T, self.Tts, self.S, self.Sts = self.load_data(isfc,mshc)

    def load_dnc(self, cdir, cfile, cvar, cnf, ckey):
        cfilel = sorted(glob.glob(cdir+'/'+self.runid+'/'+cfile))
        if len(cfilel)==0:
            print( 'no file found with this pattern '+cdir+'/'+self.runid+'/'+cfile )
            sys.exit(42)

        self.dnc['cf'+ckey]=cfilel
        self.dnc['cv'+ckey]=cvar

    def load_data(self, isfc, mshc):
        
        isf_idx=np.where(isfc.msk==1)

        # mlt distri
        area = mshc.area[isf_idx]
        drft = mshc.zisf[isf_idx]
        #
        # conversion to Gt/y
        scal = area *86400. * 365. / (1000. * 1e9) * -1.

        # melt ts
        nfile=len(self.dnc['cfmltts'])
        mltts=[] ; timmltts=[]
        _sum = np.zeros(shape=area.shape)
        for cfile in self.dnc['cfmltts']:
            timmltts.append(get_time_data(cfile,'time_centered')[:])
            mltr = get_2d_data(cfile,self.dnc['cvmltts'])[isf_idx] * scal
            mltts.append(np.sum(mltr))
            _sum = _sum + mltr

        # mlt local
        mltd = _sum/nfile 

        # mlt total
        mltt = np.sum(mltd)
      
        # sal profile
        # time series
        timprots,proSts = self.get_profilets('sal' , mshc)
        timprots,proTts = self.get_profilets('temp', mshc)

        # mean profile
        proSm=np.ma.zeros(shape=(mshc.z[:].size,2))
        proTm=np.ma.zeros(shape=(mshc.z[:].size,2))
        proSm[:,0] = mshc.z[:]
        proTm[:,0] = mshc.z[:] 
        proSm[:,1] = np.mean(proSts[:,1:],axis=1)
        proTm[:,1] = np.mean(proTts[:,1:],axis=1)

        return drft,area/(1000*1000),mltd,[timmltts, mltts],mltt,proTm,[timprots,proTts],proSm,[timprots,proSts]
    
    def get_profilets(self,cvar,mshc):
        nt = len(self.dnc['cf'+cvar+'ts'])
        nk = get_dim(self.dnc['cf'+cvar+'ts'][0],'z')
        timets=[None]*nt #np.zeros(shape=nt)
        prots=np.ma.zeros(shape=(nk,len(self.dnc['cf'+cvar+'ts'])+1))
        prots[:,0]=mshc.z[:]
        for kf,cfile in enumerate(self.dnc['cf'+cvar+'ts']) :
             prots[:,kf+1]=get_1d_data(cfile,self.dnc['cv'+cvar+'ts'])
             timets[kf]=mdates.date2num(get_time_data(cfile,'time_centered')[:][0])
        return timets,prots

    def load_att(self,runid):
        try:
            lstyle=False
            with open('STYLES/run.sty') as fid:
                for cline in fid:
                    att=cline.split('|')
                    if att[0].strip() == runid:
                        cpltname  = att[1].strip()
                        cpltclr   = att[3].strip()
                        lstyle=True
            if not lstyle:
                print( runid+' not found in STYLES/run.sty' )
                raise Exception
           
        except Exception as e:
            print( 'load_name error on ',runid )
            print( e )
            sys.exit(42)
    
        # return value
        return cpltname, cpltclr

    def __str__(self):
        return 'runid = {}, name = {}'.format(self.runid, self.name)

#================================================================================
class isf(object):
    def __init__(self, cisf, cdir='./'):
        self.name = cisf
        self.clr, self.iid, self.mapdef, self.Trange, self.Srange, self.qrange = self.load_isf(cisf)
        self.msk     = self.load_msk(cisf,cdir+'mskisf.nc','mask_isf')
        self.msk_cnt = self.load_msk(cisf,cdir+'mskisf.nc','mask_isf_front')

    def load_msk(self, cisf, cfile, cvar):
        lmsk=get_2d_data(cfile,cvar)
        lmsk[lmsk != -self.iid] = 0
        lmsk[lmsk == -self.iid] = 1
        return lmsk

    def load_isf(self, cisf):
        try:
            lstyle=False
            with open('STYLES/isf.sty') as fid:
                for cline in fid:
                    att=cline.split('|')
                    isfname      = att[0].strip()
                    if isfname == cisf :
                        clr = att[1].strip()
                        iid = int(att[2].strip())
                        mapdef = [np.asarray(eval(att[3].strip())),ccrs.Stereographic(central_latitude=-90.0, central_longitude=0)]
                        Trange= eval(att[4].strip())
                        Srange= eval(att[5].strip())
                        qrange= eval(att[6].strip())

                    lstyle=True
            if not lstyle:
                print( runid+' not found in STYLES/isf.sty' )
                raise Exception
    
            return clr,iid,mapdef,Trange,Srange,qrange

        except Exception as e:
            print( e )
            sys.exit(42)

#================================================================================
class msh(object):
    def __init__(self,cdir='./'):
        ncid   = nc.Dataset(cdir+'mask.nc')
        self.z    = ncid.variables['gdept_1d'][:].squeeze()
        ncid.close()

        self.lon  = get_2d_data(cdir+'mesh.nc','glamt')
        self.lat  = get_2d_data(cdir+'mesh.nc','gphit')
        
        self.msk  = get_2d_data(cdir+'mask.nc','tmaskutil')
        self.zisf = get_2d_data(cdir+'mask.nc','isfdraft' )
        self.area = get_2d_data(cdir+'mesh.nc','e1t'      ) * get_2d_data(cdir+'mesh.nc','e2t')

#================================================================================
class plot(object):
    def __init__(self, ctitle, cfout,isfc):
        self.ax={}
        self.title=ctitle
        self.fout=cfout
        self.fig = plt.figure(figsize=np.array([320,420])/ 25.4)

        self.ax['meltts']= self.fig.add_subplot(4, 3, (1,2))
        self.ax['mapisf']= self.fig.add_subplot(4, 3, 3, projection=isfc.mapdef[1])

        self.ax['hovT' ] = self.fig.add_subplot(4, 3, (4,5))
        self.ax['hovS' ] = self.fig.add_subplot(4, 3, (7,8))

        self.ax['TS'   ] = self.fig.add_subplot(4, 3, 9 )
        self.ax['profT'] = self.fig.add_subplot(4, 3, 10)
        self.ax['profS'] = self.fig.add_subplot(4, 3, 11)

        self.ax['aread'] = self.fig.add_subplot(4, 3, 12)
        self.ax['meltd'] = self.ax['aread'].twiny()

        plt.subplots_adjust(left=0.1, right=0.98, bottom=0.08, top=0.92, wspace=0.15, hspace=0.26)

    def plot_dist(self,runlst,isfc,obsc,mshc):
        
        # melt distribution
        self.plot_aread(runlst)
        self.plot_meltd(runlst)

        # total melt ts
        self.plot_meltts(runlst,obsc)

        # obs temperature profile
        self.plot_obsTS(obsc)
        self.plot_TS(runlst)
       
        self.plot_obs('profT',obsc.name,obsc.obsT,obsc.Tz)
        self.plot_obs('profS',obsc.name,obsc.obsS,obsc.Sz)

        for irun, runid in enumerate(runlst):
            self.plot_prof('profT',runid,runid.Tts,runid.T)
            self.plot_prof('profS',runid,runid.Sts,runid.S)

        if len(runlst) == 1:
            # temperture hovmuller
            self.plot_hov('hovT',runlst[0].Tts,isfc.Trange[0],isfc.Trange[1],'[C]')
            self.plot_hov('hovS',runlst[0].Sts,isfc.Srange[0],isfc.Srange[1],'[g/kg]')

        # map isf
        self.plot_mapisf(isfc,obsc,mshc)

        self.add_legend()

        self.add_title()

        self.set_axes('meltts', 'bottom',                                         cylabel='Total melt [Gt/y]', ylim=isfc.qrange, lyrevert=False)
        self.set_axes('hovT'  , 'bottom', lxticklabels=False,                     cylabel='Depth [m]'        , ylim=[0, 1000])
        self.set_axes('hovS'  , 'bottom',                                         cylabel='Depth [m]'        , ylim=[0, 1000])
        self.set_axes('profT' , 'bottom',                                         cylabel='Depth [m]'        , ylim=[0, 1000]  , xlim=isfc.Trange)
        self.set_axes('profS' , 'bottom',                     lyticklabels=False,                              ylim=[0, 1000]  , xlim=isfc.Srange)
        self.set_axes('TS'    , 'bottom',                                                                      ylim=isfc.Trange, xlim=isfc.Srange, lyrevert=False) 
        self.set_axes('aread' , 'top'   , lxticklabels=False, lyticklabels=False,                              ylim=[0, 1000]  , lxtick=False)
        self.set_axes('meltd' , 'bottom',                     lyticklabels=False,                              ylim=[0, 1000])
        self.set_axes('mapisf', 'bottom', lxticklabels=False,                                                                  lyrevert=False)

        self.save()

    def plot_mapisf(self,isfc,obsc,mshc):
        isf_features   = cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m',facecolor='none',edgecolor='k')
        coast_features = cartopy.feature.NaturalEarthFeature('physical', 'coastline'                      , '50m',facecolor='0.75',edgecolor='k')

        cmap=plt.get_cmap('Blues',4)
        cmap.set_under('w', 1.0)
        coords = isfc.mapdef[1].transform_points( ccrs.PlateCarree(), isfc.mapdef[0][0:2], isfc.mapdef[0][2:4])
        self.ax['mapisf'].set_extent([coords[0, 0], coords[1, 0], coords[0, 1], coords[1, 1]], isfc.mapdef[1])
        self.ax['mapisf'].gridlines()
        for feature in [coast_features, isf_features]:#, bathy1000_features, bathy2000_features, bathy3000_features]:
            self.ax['mapisf'].add_feature(feature)
        mshc.lon[mshc.lon<0]=mshc.lon[mshc.lon<0]+360.
        self.ax['mapisf'].pcolormesh(mshc.lon[0:-4,:],mshc.lat[0:-4,:],isfc.msk[0:-4,:],vmin=0.5,vmax=2,cmap=cmap,transform=ccrs.PlateCarree(),rasterized=True, label=isfc.name, shading='flat')
        self.ax['mapisf'].plot(obsc.lon,obsc.lat,'o', markersize=6, color='royalblue', markeredgecolor='royalblue', transform=ccrs.PlateCarree(), label='WOA2018 obs')

    def plot_hov(self,cax,dat,rmin,rmax,cunit):
        # plot
        pcol=self.ax[cax].pcolormesh(dat[0][:],dat[1][:,0],dat[1][:,1::],vmin=rmin, vmax=rmax, shading='nearest')
        #
        # xaxis
        datemin = mdates.num2date(dat[0][0])
        datemax = mdates.num2date(dat[0][-1])
        nyear=(datemax-datemin).days/365 + 1
        if nyear < 10:
            nyt=1
        elif 10<=nyear<50:
            nyt=5
        elif 50<=nyear<100:
            nyt=10
        else:
            nyt=100
        self.ax[cax].xaxis.set_major_locator(mdates.YearLocator(nyt,month=1,day=1))
        self.ax[cax].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        #
        # cbar
        cbar = plt.colorbar(pcol, ax=self.ax[cax], extend='both')
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_title(cunit,fontsize=16)

    def plot_obsTS(self,obsc):
        TScastid = set(obsc.Scast).intersection(obsc.Tcast)
        idxT=[list(obsc.Tcast).index(iid) for iid in TScastid]
        idxS=[list(obsc.Scast).index(iid) for iid in TScastid]
        for icast in range(len(idxT)):
            self.ax['TS'].plot(obsc.obsS[:,idxS[icast]],obsc.obsT[:,idxT[icast]],color='gainsboro',linewidth=2,alpha=0.6) 
        self.ax['TS'].plot(np.nanmean(obsc.obsS,axis=1),np.nanmean(obsc.obsT,axis=1),color='silver',linewidth=5,label=obsc.name,zorder=99)

    def plot_TS(self,runlst):
        for irun, runid in enumerate(runlst):
            for jp in range(0,len(runid.Tts[0])):
                self.ax['TS'].plot(runid.Sts[1][:,jp+1],runid.Tts[1][:,jp+1],color=lighten_color(runid.clr,0.5),linewidth=2,alpha=0.6)
            self.ax['TS'].plot(runid.S[:,1],runid.T[:,1],color=runid.clr,linewidth=5,label=runid.name,zorder=99)

        # ======================= TS diag ====================================
        ydim=xdim=100
        dens = np.zeros((ydim,xdim))
        ti = np.linspace(-3.0,5.0 ,ydim)
        si = np.linspace(33.4,35.4,xdim)
        for j in range(0,int(ydim)):
            for i in range(0, int(xdim)):
                dens[j,i]=gsw.density.sigma0(si[i], ti[j])
        cnt_lvl=np.arange(26.0,29.0,0.2)
        CS = self.ax['TS'].contour(si, ti, dens, cnt_lvl, linestyles='dashed', colors='0.7',zorder=99, linewidths=1)
        self.ax['TS'].clabel(CS, inline=True, fontsize=12, fmt='%4.2f', inline_spacing=1) # Label every second level

    def plot_prof(self,cax,runid,datts,dat):
        datmin=np.amin(datts[1][:,1:],axis=1)
        datmax=np.amax(datts[1][:,1:],axis=1)
        self.ax[cax].fill_betweenx(datts[1][:,0],datmin[:],datmax[:],color=lighten_color(runid.clr,0.5),alpha=0.5)
        self.ax[cax].plot(dat[:,1],dat[:,0],color=runid.clr,linewidth=5,label=runid.name)

    def plot_obs(self,cax,clabel,obsdat,obsz):
        datmin=np.nanmin(obsdat[:,:],axis=1)
        datmax=np.nanmax(obsdat[:,:],axis=1)
        zdat  =np.nanmean(obsz[:,:],axis=1)
        self.ax[cax].fill_betweenx(zdat,datmin,datmax,color='gainsboro',alpha=0.5)
        self.ax[cax].plot(np.nanmean(obsdat,axis=1),np.nanmean(obsz,axis=1),color='silver',linewidth=5,label=clabel)

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
        timemax=np.max([runlst[irun].meltts[0][-1] for irun in range(0,len(runlst))])
        timemin=np.min([runlst[irun].meltts[0][0] for irun in range(0,len(runlst))])
        timemax=timemax + (timemax-timemin)/20.

        for irun, runid in enumerate(runlst):
            self.ax['meltts'].plot(runid.meltts[0][:], runid.meltts[1][:], '-o', markeredgecolor=runid.clr, color=runid.clr, markersize=8)
        self.ax['meltts'].errorbar(timemax, obsc.obsm[0], yerr=obsc.obsm[1], fmt='o', markersize=8, color='k', linewidth=3, label='Rignot et al. 2013')

    def add_title(self):
        self.ax['hovT'  ].set_title('b) Temperature'         ,fontsize=18)
        self.ax['hovS'  ].set_title('c) Salinity'            ,fontsize=18)
        self.ax['TS'    ].set_title('d) TS diagram'          ,fontsize=18)
        self.ax['profT' ].set_title('e) Temperature [C]'     ,fontsize=18)
        self.ax['profS' ].set_title('f) Salinity [g/kg]'     ,fontsize=18)
        self.ax['meltd' ].set_title('g) Ice shelf melt [G/y]',fontsize=18)
        self.ax['meltts'].set_title('a) Ice shelf melt [G/y]',fontsize=18)
        self.fig.suptitle(self.title,fontsize=24)

    # set axes bound and grid
    def set_axes(self, cname, cxpos, lxticklabels=True, lyticklabels=True, cxlabel=None, cylabel=None, lyrevert=True, xlim=None, ylim=None, lxtick=True):
    #def set_axes(self, cname, lxtick, lxticklabel, cxlabel, cxpos, cylabel=None, lyrevert=True, ylim=None, lxtick=True):
        if ylim:
            self.ax[cname].set_ylim(ylim)

        if xlim:
            self.ax[cname].set_xlim(xlim)

        if lyrevert:
            self.ax[cname].invert_yaxis()

        self.ax[cname].tick_params(labelsize=16)
        self.ax[cname].xaxis.set_ticks_position(cxpos)
        self.ax[cname].grid(ls=':',c='k')
        self.ax[cname].set_axisbelow(False)

        if not lxtick:
            self.ax[cname].set_xticks([])
        else:
            if not lxticklabels:
                self.ax[cname].set_xticklabels([])

        if cylabel:
            self.ax[cname].set_ylabel(cylabel, fontsize=16)

        if not lyticklabels:
            self.ax[cname].set_yticklabels([])

        if cxlabel:
            self.ax[cname].set_xlabel(cxlabel, fontsize=16)

    # save plot
    def save(self):
        self.fig.savefig('FIGURES/'+self.fout+'.png', format='png', dpi=150)

    def add_legend(self, ncol=5, lvis=True):
        # isf area
        lline, llabel = self.ax['aread'].get_legend_handles_labels()
        leg = self.ax['aread'].legend(lline, llabel, loc='lower right', fontsize=16, frameon=False)
        for item in leg.legendHandles:
            item.set_visible(lvis)

        # profile obs
        lline, llabel = self.ax['profT'].get_legend_handles_labels()
        leg=  self.ax['profT'].legend([lline[0]], [llabel[0]], loc='lower right', fontsize=16, frameon=False)

        # isf melt
        lline, llabel = self.ax['meltts'].get_legend_handles_labels()
        leg=  self.ax['meltts'].legend([lline[-1]], [llabel[-1]], numpoints=1, loc='lower left', fontsize=16, frameon=False)

        # model profile
        lline, llabel = self.ax['profT'].get_legend_handles_labels()
        lax = plt.axes([0.00, 0.00, 1, 0.06])
        leg=  plt.legend(lline[1::], llabel[1::], loc='upper left', ncol = ncol, fontsize=16, frameon=False)
        for item in leg.legendHandles:
            item.set_visible(lvis)
        lax.set_axis_off()

# ================================= NETCDF ===============================
def get_name(regex,varlst):
    revar = re.compile(r'\b%s\b'%regex,re.I)
    cvar  = list(filter(revar.match, varlst))
    if (len(cvar) > 1):
        print( regex+' name list is longer than 1 or 0; error' )
        print( cvar )
        print( cvar[0]+' is selected' )
    if (len(cvar) == 0):
        print( 'no match between '+regex+' and :' )
        print( varlst )
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
        print( 'cshape undefined, error' )
        print( dimlst )
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
        print( 'dimension direction unknown, need to be x, y, z or k' )
        sys.exit(42)

    cdim=list(filter(redim.match, ncid.dimensions.keys()));
    if (len(cdim) > 1):
        print( regex+' name list is longer than 1; error' )
        print( cdim )
        sys.exit(42)
    elif (len(cdim) == 0):
        #print( cdir+' dim in '+cfile+' is 0.' )
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
            timeidx[itime] = np.datetime64(time[itime])#,'us')
        else:
            timeidx[itime] = np.datetime64(time)       #,'us')

    return timeidx[:]

# get_1d_data
def get_1d_data(cfile,cvar,ktime=0,kx=0,ky=0):
    print( ' reading '+cvar+' in '+cfile+' ...' )
 
    nx,ny,nz,nt=get_dims(cfile)

    if ((nx > 1) or (ny > 1) or (nt > 1)):
        print('data format not supported, dimension along t or x or y larger than 1')
        sys.exit(42)

    ncid  = nc.Dataset(cfile)
    clvar = get_name(cvar,ncid.variables.keys())
    var   = ncid.variables[clvar]
    shape = get_variable_shape(ncid,var)

    dat1d=var[:].squeeze()    

    ncid.close()
    return dat1d 

# get_2d_data
def get_2d_data(cfile,cvar,ktime=0,klvl=0):
    print( ' reading '+cvar+' in '+cfile+' ...' )

    nx,ny,nz,nt=get_dims(cfile)

    ncid  = nc.Dataset(cfile)
    clvar = get_name(cvar,ncid.variables.keys())
    var   = ncid.variables[clvar]
    shape = get_variable_shape(ncid,var)

    if shape=='XY' :
        if (klvl > 0) :
            print( 'error klvl larger than 0 (klvl = '+str(klvl)+')' )
            sys.exit(42)
        if (ktime > 0) :
            print( 'error ktime larger than 0 (ktime = '+str(ktime)+')' )
            sys.exit(42)
        dat2d=var[:,:]
    elif shape=='XYT' :
        if (klvl > 0) :
            print( 'error klvl larger than 0 (klvl = '+str(klvl)+')' )
            sys.exit(42)
        if (ktime > 0) :
            print( 'error ktime larger than 0 (ktime = '+str(ktime)+')' )
            sys.exit(42)
        dat2d=var[ktime,:,:]
    elif shape=='XYZ' :
        if (ktime > 0) :
            print( 'error ktime larger than 0 (ktime = '+str(ktime)+')' )
            sys.exit(42)
        dat2d=var[klvl,:,:]
    elif len(var.shape)==4 :
        print( ' 4d variable XYZT' )
        dat2d=var[ktime,klvl,:,:]
    else:
        print( cvar+' contains '+str(len(var.shape))+' dimensions' )
        print( 'dimension names are '+var.dimensions )
        print( ' shape unknown, exit ' )
        sys.exit(42)
    ncid.close()
    return dat2d

#================================================================================
# check argument
def load_argument():
    # deals with argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-runid", metavar='runid list'   , help="used to look information in runid.db"                  , type=str, nargs="+" , required=True )

    parser.add_argument("-fmltts", metavar='file list ts', help="file list to plot (default is runid_var.nc)"           , type=str, nargs="+" , required=False)
    parser.add_argument("-ftemts" , metavar='file T list'  , help="file list to plot (default is runid_var.nc)"           , type=str, nargs="+" , required=False)
    parser.add_argument("-fsalts" , metavar='file T list'  , help="file list to plot (default is runid_var.nc)"           , type=str, nargs="+" , required=False)

    parser.add_argument("-vmlt" , metavar='isf mlt var'  , help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs="+" , required=True )
    parser.add_argument("-vtem" , metavar='oce sal var'  , help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs="+" , required=True )
    parser.add_argument("-vsal" , metavar='oce sal var'  , help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs="+" , required=True )

    parser.add_argument("-isf"  , metavar='isf name'     , help="variable to look for in the netcdf file ./runid_var.nc", type=str, nargs=1   , required=True )

    parser.add_argument("-dir"  , metavar='directory of input file' , help="directory of input file"                    , type=str, nargs=1   , required=False, default=['./'])
    parser.add_argument("-title", metavar='figure title', help="figure title"                                           , type=str, nargs=1   , required=True)
    parser.add_argument("-o"    , metavar='figure_name', help="output figure name without extension"                    , type=str, nargs=1   , required=False, default=['output'])
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
       print( 'all argument ('+ctxt+') set to the same value' )
    elif (len(argin) == nrun):
       argout = argin
    else:
       print( 'number of arguments is not compatible with the number of runid ('+ctxt+')' )
       sys.exit(42)
    return argout
#==============================================================================

def main():

# load argument
    args = load_argument()

# output argument list
    output_argument_lst('FIGURES/'+args.o[0]+'.txt', sys.argv)

# initialisation
    print('initialisation')
    nrun=len(args.runid)
    run_lst=[None]*nrun

#    fmlt  =fix_arglist(nrun,args.fmlt  ,'args.fmlt')
    vmlt  =fix_arglist(nrun,args.vmlt  ,'args.vmlt')
#    ftem  =fix_arglist(nrun,args.ftem  ,'args.ftem')
    vtem  =fix_arglist(nrun,args.vtem  ,'args.vtem')
#    fsal  =fix_arglist(nrun,args.fsal  ,'args.fsal')
    vsal  =fix_arglist(nrun,args.vsal  ,'args.vsal')
    fmltts=fix_arglist(nrun,args.fmltts,'args.fmltts')
    ftemts=fix_arglist(nrun,args.ftemts,'args.ftemts')
    fsalts=fix_arglist(nrun,args.fsalts,'args.fsalts')
    cdir  =fix_arglist(nrun,args.dir   ,'args.dir')

# load data for each run
    print('load each run')
    for irun, runid in enumerate(args.runid):

        # load isf data
        print('load isf data')
        isfc = isf(args.isf[0],cdir[irun]+'/'+runid+'/')

        # load msh data
        print('load msh')
        mshc = msh(cdir[irun]+'/'+runid+'/')

        # load WOA profile
        print('load profile')
        profc = obs(isfc,mshc)

        #run_lst[irun] = run(runid, cdir[irun], ftem[irun], fsal[irun], fmlt[irun], ftemts[irun], fsalts[irun], fmltts[irun], vtem[irun], vsal[irun], vmlt[irun], isfc, mshc)
        run_lst[irun] = run(runid, cdir[irun], ftemts[irun], fsalts[irun], fmltts[irun], vtem[irun], vsal[irun], vmlt[irun], isfc, mshc)

    print('plot')
    plt_isf = plot(args.title[0], args.o[0],isfc)
    plt_isf.plot_dist(run_lst,isfc,profc,mshc)

    if args.noshow:
       pass
    else:
       plt.show()

if __name__=="__main__":
    main()
