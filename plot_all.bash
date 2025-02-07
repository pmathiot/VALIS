#echo 'Plot time series'
#python2.7 plot_mlt_timeseries.py -dir $DIR -runid $RUNID -f isfmlt_ts_u-*.nc -var TOTA -title 'ISF melt (all isf)' -o melt_ts_total -obs Rignot_2013 -minmax 0 3500 -sf -1 -noshow
#python2.7 plot_mlt_timeseries.py -dir $DIR -runid $RUNID -f isfmlt_ts_u-*.nc -var AMER ROSS FRIS LARC FIMB -title 'ISF melt (cold isf)' -o melt_ts_cold -obs Rignot_2013 -minmax 0 600 -sf -1 -noshow
#python2.7 plot_mlt_timeseries.py -dir $DIR -runid $RUNID -f isfmlt_ts_u-*.nc -var SLUZ GETZ PINE TWAI GEVI -title 'ISF melt (warm isf)' -o melt_ts_warm -obs Rignot_2013 -minmax 0 1000 -sf -1 -noshow
#wait

#echo 'Plot map'
#RUNID=bt705
#FILET="nemo_${RUNID}o_1y_${TAG}_grid-T.nc"
#python2.7 plot_meltrate_sector.py -ftem ${DIR}/${RUNID}/${FILET} -fisf ${DIR}/${RUNID}/${FILET} -vtem thetaob -visf sowflisf -t "GO8 JRA (${YEAR})" -o melt_sector_${RUNID}_${YEAR}.png

#RUNID=bt036
#FILET="nemo_${RUNID}o_1y_${TAG}_grid-T.nc"
#python2.7 plot_meltrate_sector.py -ftem ${DIR}/${RUNID}/${FILET} -fisf ${DIR}/${RUNID}/${FILET} -vtem thetaob -visf sowflisf -t "GO8 DFS5.2 (${YEAR})" -o melt_sector_${RUNID}_${YEAR}.png

#RUNID=bt223
#FILET="nemo_${RUNID}o_1y_${TAG}_grid-T.nc"
#python2.7 plot_meltrate_sector.py -ftem ${DIR}/${RUNID}/${FILET} -fisf ${DIR}/${RUNID}/${FILET} -vtem thetaob -visf sowflisf -t "GO8 CORE (${YEAR})" -o melt_sector_${RUNID}_${YEAR}.png

if [ $# == 0 ] ; then 'usage: plot_all.bash [RUNID (CONFIG-CASE)] ; ecit 42'; exit 42 ; fi

RUNID=$1
FILET_TS="ISF_Tprof*1y_y*_gridT.nc"
FILES_TS="ISF_Sprof*1y_y*_gridT.nc"
FILEQ_TS="*y????.1y*flxT.nc"

. param.bash

if [ $RUN_MLT_DST == 1 ] ; then 
   for ISF in FIMB RIIS PINE FRIS ROSS AMER GETZ GEVI ; do
       echo "$ISF ..."
       python SCRIPTS/plot_mlt_distri.py -runid ${RUNID} -fmltts ${FILEQ_TS} -vmlt sowflisf_cav -ftemts ${FILET_TS} -vtem mean_votemper_${ISF} -fsalts ${FILES_TS} -vsal mean_vosaline_${ISF} -isf ${ISF} -dir ${WRKPATH} -title "${ISF}" -o ${ISF} -noshow
   done
else
   echo '$RUN_MLT_DST is different to 1; nothing is done'
fi
wait
