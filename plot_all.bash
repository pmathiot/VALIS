RUNID='bt036 bt223 bt705'
TAG=1996-2005 ; YEAR=1996-2005
DIR=/scratch/pmathiot/MLT

echo 'Plot time series'
python2.7 plot_mlt_timeseries.py -dir $DIR -runid $RUNID -f isfmlt_ts_u-*.nc -var TOTA -title 'ISF melt (all isf)' -o melt_ts_total -obs Rignot_2013 -minmax 0 3500 -sf -1 -noshow
python2.7 plot_mlt_timeseries.py -dir $DIR -runid $RUNID -f isfmlt_ts_u-*.nc -var AMER ROSS FRIS LARC FIMB -title 'ISF melt (cold isf)' -o melt_ts_cold -obs Rignot_2013 -minmax 0 600 -sf -1 -noshow
python2.7 plot_mlt_timeseries.py -dir $DIR -runid $RUNID -f isfmlt_ts_u-*.nc -var SLUZ GETZ PINE TWAI GEVI -title 'ISF melt (warm isf)' -o melt_ts_warm -obs Rignot_2013 -minmax 0 1000 -sf -1 -noshow
wait

echo 'Plot map'
RUNID=bt705
FILET="nemo_${RUNID}o_1y_${TAG}_grid-T.nc"
python2.7 plot_meltrate_sector.py -ftem ${DIR}/${RUNID}/${FILET} -fisf ${DIR}/${RUNID}/${FILET} -vtem thetaob -visf sowflisf -t "GO8 JRA (${YEAR})" -o melt_sector_${RUNID}_${YEAR}.png

RUNID=bt036
FILET="nemo_${RUNID}o_1y_${TAG}_grid-T.nc"
python2.7 plot_meltrate_sector.py -ftem ${DIR}/${RUNID}/${FILET} -fisf ${DIR}/${RUNID}/${FILET} -vtem thetaob -visf sowflisf -t "GO8 DFS5.2 (${YEAR})" -o melt_sector_${RUNID}_${YEAR}.png

RUNID=bt223
FILET="nemo_${RUNID}o_1y_${TAG}_grid-T.nc"
python2.7 plot_meltrate_sector.py -ftem ${DIR}/${RUNID}/${FILET} -fisf ${DIR}/${RUNID}/${FILET} -vtem thetaob -visf sowflisf -t "GO8 CORE (${YEAR})" -o melt_sector_${RUNID}_${YEAR}.png

echo 'Plot isf sanity plot'
RUNID='bt036 bt223 bt705 bs852'
FILES="nemo_*o_1y_${TAG}_grid-T.nc"
FILETS="nemo_*o_1y_*01-*01_grid-T.nc"
VARMLT='sowflisf'
VART='thetao'
for ISF in PINE RIIS FIMB FRIS ROSS AMER GETZ GEVI ; do #FIMB RIIS LARC AMER ROSS PINE GEVI GETZ; do
   TITLE="$ISF ($YEAR)"
   python2.7 plot_mlt_distri.py -dir $DIR -runid $RUNID -fmlt $FILES -fmltts $FILETS -vmlt $VARMLT -vtem $VART -isf $ISF -title "$TITLE" -o ${ISF}_${YEAR} -noshow > log_$ISF &
done
wait

