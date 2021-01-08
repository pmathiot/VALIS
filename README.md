# VALIS

## Purpose
<!--
 * Assess the order 0 of the water masses properties and ice shelf melt in NEMO output
   * plot_meltrate_sector.py plots bottom temperature and ice shelf melt map per sector for a specific run
![Alt text](melt_sector.png?raw=true "Example of plot_meltrate_sector.py output")

   * plot_mlt_timeseries.py plots ice shelf melt time series against observation/model estimates
![Alt text](melt_ts_cold.png?raw=true "Example of plot_mlt_timeseries.py output")
-->
   * **plot_mlt_distri.py** assess for a specific ice shelf :
      * water mass properties at the calving front against WOA data
      * ice shelf melt distribution vs isf draft area distribution
      * ice shelf melt time serie
      * T and S hovmuller plot at calving front
      * comparison between run profile and WOA2018 profile
      * localisation map
![Alt text](FIGURES/FRIS.png?raw=true "Example of plot_mlt_distri.py output")

## Usage

* Add your WOA2018 profile in OBS:
   - file name 'ISF name'_CTD_WOA2018.nc. 
   - You can generate this file by going here: https://www.nodc.noaa.gov/OC5/SELECT/dbsearch/dbsearch.html

* set up your style in STYLE:
   - define your color, line style, name (...) for every simulation        in run.sty
   - define your color, line style, id   (...) for every ice shelves       in isf.sty 
   - define your color, line style, name (...) for every melt rate sources in obs.sty

* if required, add your melt rate comparison data in MLT

* setup the param file:
   - WRKPATH need to be the same as the one used in VALSO (ie the master directory where all the simulation sub directories are)
   - define the kind of plot you want

* run ./plot_all.bash [RUNID]
   - RUNID have a syntax like this for example (eORCA025.L121-OPM006)
   - all available data are processed by default

<!--
 * plot_mlt_timeseries.py example: ```python2.7 plot_mlt_timeseries.py -dir [DATA DIR] -runid [RUNID/NAME] -f [FILENAMES (wildcard accepted)] -var [ISF list] -title [TITLE] -o [OUTPUT name] -obs [OBS file] -minmax [DATA range] -sf [SCALE factor] -noshow```
    * DATA DIR: where the data are (expect an dir tree like this DATADIR/RUNID/*.nc)
    * RUNID: simulation key work used to retreive the run style in run.sty and the data directory (1 arg per run)
    * FILENAMES: list of all the files needed (wild card accepted). Each file should have been computed using cdfisf_diags (see CDFTOOLS repository) before.
    * ISF list: list of ice shelf name (see in isf.sty for accepted name)
    * TITLE: plot title
    * OUTPUT name: output name without the extension (a png and txt file will be saved). The txt file contain the exact python command run and the png the figure.
    * OBS: obs file name (see Rignot_2013.txt for template)
    * DATA range: figure y range 
    * SCALE factor: scale factor to apply to the data
    * noshow: only save the output (no display)

 * plot_meltrate_sector.py example: ```python2.7 plot_meltrate_sector.py -ftem [FILET] -fisf [FILEISF] -vtem [BOTTOMT var] -visf [ISFMLT var (kg/m2/s)] -t [TITLE] -o [OUTPUT]```
    * FILET  : netcdf file containing [BOTTOMT var] variable
    * FILEISF: netcdf file containing [FILEISF] (ice shelf melt rate) variable
    * TITLE: plot title
    * OUTPUT: output figure name
-->
## Details

 ### RUN_MLT_DST: plot_mlt_distri.py

```
   python plot_mlt_distri.py -dir [DATA DIR] -runid [RUNID] -fmltts [FILEMLT] -vmlt [VARISF] -ftemts [FILET] -vtem [VART] -fsalts [FILES] -vsal [VARS] -isf [ISF] -title "$TITLE" -o ${ISF}_${YEAR} -noshow
```

 * `DATADIR`: where the data are (expect a dir tree like this DATADIR/RUNID/\*.nc)
 * `RUNID`  : simulation key work used to retreive the run style in run.sty and the data directory (1 arg per run)
 * `FILEMLT`: netcdf file containing integrated melt per isf.
 * `FILET`  : netcdf file containing mean T profile in front of each ice shelf
 * `FILES`  : netcdf file containing mean S profile in front of each ice shelf
 * `ISF`    : ice shelf name (see isf.sty for accepted name)
 * `TITLE`  : plot title
 * `OUTPUT` : output file name without the extension (a png and txt file will be saved). The txt file contain the exact python command run and the png the figure.
 * `noshow` : only save the output (no display)
    
## Requirements:
 * python with the following module:
```
name: py37
channels:
  - defaults
dependencies:
  - python=3.7.4
  - cartopy
  - gsw
  - scipy
  - netcdf4
  - dask
  - xarray
prefix: /home/pmathiot/.conda/envs/py37
```
 * plot_mlt_distri.py and plot_meltrate_sector.py requires to have all the run to use the same mesh and mask files.
 * obs.sty file: you need to add the obs file you want to plot in the list of available obs style.

## Recommandations:
 * Run VALSO toolbox before (at least the ISF diagnostics) to have all the input file needed for the plot processed with the right format and everything ...
