import matplotlib.pyplot as plt
import scipy.io
from numpy import ma
import matplotlib as matplotlib
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy
import numpy as np
from netCDF4 import Dataset
from matplotlib import gridspec
matplotlib.use('GTKAgg')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-ftem", metavar='file_name_T2d' , help="input file name for temp2d"  , type=str  , nargs=1, required=False)
parser.add_argument("-fisf", metavar='file_name_melt', help="input file name for isf melt", type=str  , nargs=1, required=True)
parser.add_argument("-vtem", metavar='var_name_T2d'  , help="var name for temp2d"         , type=str  , nargs=1, required=False)
parser.add_argument("-visf", metavar='var_name_melt' , help="var name for isf melt"       , type=str  , nargs=1, required=True)
parser.add_argument("-o"   , metavar='output_file'   , help="output figure name (png)"    , type=str  , nargs=1, required=True)
parser.add_argument("-t"   , metavar='figure_title'  , help="figure title"                , type=str  , nargs=1, required=True)
args = parser.parse_args()

test=scipy.io.loadmat('/data/users/pmathiot/BCKUP/BAS_BCKUP/MATLAB/PLOT_PERIANT025/meltratepal.mat')
rgb=test["meltratepal"]
cmap=mcolors.ListedColormap(rgb)
cmap.set_bad('0.75', 1.0)
cmap_tem=plt.get_cmap('RdBu_r',20)

plt.figure(figsize=np.array([200,300]) / 25.4)

imin=1 ; imax=1440
jmin=0 ; jmax=600

ncid   = Dataset('coordinates.nc')
lat2d  = ncid.variables['gphit'][0,jmin:jmax,imin:imax]
lon2d  = ncid.variables['glamt'][0,jmin:jmax,imin:imax]
ncid.close()

ncid   = Dataset('mesh_mask.nc')
msk    = ncid.variables['tmaskutil'][0,jmin:jmax,imin:imax]
msk_surf = ncid.variables['tmask'][0,0,jmin:jmax,imin:imax]
ncid.close()

# isf melt
ncid   = Dataset(args.fisf[0])
var2d  = ncid.variables[args.visf[0]][0,jmin:jmax,imin:imax]
var2d  = -1*var2d * 86400 * 365/1000
var2dm = ma.masked_where(msk*var2d==0.0,var2d)
ncid.close()

# ocean temperature
ltemp=False
if args.ftem:
    ncid   = Dataset(args.ftem[0])
    temp2d  = ncid.variables[args.vtem[0]][0,jmin:jmax,imin:imax]
    temp2dm = ma.masked_where(msk_surf*temp2d==0.0,temp2d)
    ncid.close()
    ltemp=True

isf_features   = cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines', '50m',facecolor='none',edgecolor='k')
coast_features = cartopy.feature.NaturalEarthFeature('physical', 'coastline'                         , '50m',facecolor='0.75',edgecolor='k')

ax = [None] * 4
bc = [None] * 4

gs = gridspec.GridSpec(3, 2)
plot_proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=-100.0)
ax[0] = plt.subplot(gs[0,:], projection=plot_proj)
ax[0].set_extent([-1.354e6, 1.301e6, 1.030e6, 2.215e6],plot_proj)
ax[0].add_feature(coast_features)
ax[0].add_feature(isf_features)
ax[0].gridlines()
print lon2d.shape, lat2d.shape, var2dm.shape
pcol_hmlt = ax[0].pcolormesh(lon2d,lat2d,var2dm ,vmin=-6.5 ,vmax=30,cmap=cmap    ,transform=ccrs.PlateCarree(),rasterized=True)
if ltemp:
    pcol_tem  = ax[0].pcolormesh(lon2d,lat2d,temp2dm,vmin=-2   ,vmax=2 ,cmap=cmap_tem,transform=ccrs.PlateCarree(),rasterized=True)
ax[0].set_title('(a) West Antarctica ice shelves',fontsize=16)

plot_proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=-180.0)
ax[1] = plt.subplot(gs[1,0], projection=plot_proj)
ax[1].set_extent([-5.115e5, 6.499e5, 4.591e5, 1.6205e6],plot_proj)
ax[1].add_feature(coast_features)
ax[1].add_feature(isf_features)
ax[1].gridlines()
ax[1].pcolormesh(lon2d,lat2d,var2dm ,vmin=-1.1 ,vmax=5 ,cmap=cmap    ,transform=ccrs.PlateCarree(),rasterized=True)
if ltemp:
    ax[1].pcolormesh(lon2d,lat2d,temp2dm,vmin=-2   ,vmax=2 ,cmap=cmap_tem,transform=ccrs.PlateCarree(),rasterized=True)
print np.max(msk[:,1])
ax[1].set_title('(b) Ross ice shelf',fontsize=16)

plot_proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=-30.0)
ax[2] = plt.subplot(gs[1,1], projection=plot_proj)
ax[2].set_extent([-1.507e6, 6.97e5, 5.872e5, 2.7912e6],plot_proj)
ax[2].add_feature(coast_features)
ax[2].add_feature(isf_features)
ax[2].gridlines()
ax[2].pcolormesh(lon2d,lat2d,var2dm ,vmin=-1.1 ,vmax=5 ,cmap=cmap    ,transform=ccrs.PlateCarree(),rasterized=True)
if ltemp:
    ax[2].pcolormesh(lon2d,lat2d,temp2dm,vmin=-2   ,vmax=2 ,cmap=cmap_tem,transform=ccrs.PlateCarree(),rasterized=True)
print np.max(msk[:,1])
ax[2].set_title('(c) Weddell Sea ice shelves',fontsize=16)

plot_proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=-90.0)
ax[3] = plt.subplot(gs[2,:], projection=plot_proj)
ax[3].set_extent([-2.480e6, 2.526e6, -2.944e6, -6.308e5],plot_proj)
ax[3].coastlines(resolution='50m')
ax[3].add_feature(coast_features)
ax[3].add_feature(isf_features)
ax[3].gridlines()
pcol_lmlt = ax[3].pcolormesh(lon2d,lat2d,var2dm ,vmin=-1.1 ,vmax=5 ,cmap=cmap    ,transform=ccrs.PlateCarree(),rasterized=True)
if ltemp:
    pcol_tem  = ax[3].pcolormesh(lon2d,lat2d,temp2dm,vmin=-2   ,vmax=2 ,cmap=cmap_tem,transform=ccrs.PlateCarree(),rasterized=True)
print np.max(msk[:,1])
ax[3].set_title('(d) East Antarctica ice shelves',fontsize=16)

plt.subplots_adjust(left=0.01, right=0.80, bottom=0.01, top=0.95, wspace=0.1, hspace=0.1)
for iplt in range(0,4):
    ax[iplt].apply_aspect()
    bc[iplt] = ax[iplt].get_position()

cax = plt.axes([bc[0].x1+0.02, bc[0].y0, 0.02, bc[0].y1-bc[0].y0])
vlevel=[-5, 0, 5, 10, 15, 20, 25, 30];
cbar= plt.colorbar(pcol_hmlt, ticks=vlevel, cax=cax)
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_title('m/y',fontsize=14)

cax = plt.axes([bc[2].x1+0.02, bc[3].y0, 0.02, bc[2].y1-bc[3].y0])
vlevel=[-1, 0, 1, 2, 3, 4, 5];
cbar= plt.colorbar(pcol_lmlt, ticks=vlevel, cax=cax)
cbar.ax.tick_params(labelsize=14)

if ltemp:
    cax = plt.axes([bc[3].x1+0.1, bc[3].y0, 0.02, bc[0].y1-bc[3].y0])
    cbar= plt.colorbar(pcol_tem, cax=cax)
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_title('$\degree$C',fontsize=14)

plt.suptitle(args.t[0],fontsize=16)

plt.savefig(args.o[0], format='png', dpi=300)

plt.show()
