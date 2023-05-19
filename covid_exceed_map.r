# Script Title: covid_exceed_map.r
# Author: Massimo Cavallaro
# Description: Plot color maps for each week recursively from the warning
# scores dumped by the script `rancover_exceed.r`.
# Dependencies: read_data.r, create_matrices.r

source('covid_exceed_map_utils.r')

# day_names = as.POSIXct('2022-2-21', '%Y-%m-%d', tz='GMT') - as.difftime(782:0,units='days')
  
source("ws_250-783_2022-08-31_mean_rincidence.R")

for (i in 21:532){
  figure_name = sprintf("Exceed_map_%s.png", sprintf("%03d", i))
  png(figure_name,
    width = 4 * 3,
    height = 6,
    pointsize = 12,
    unit = 'in',
    res = 600)
  plot_exceed_map3(
    ws.p = ws.p,
    ws.f = ws.f,
    averages.p = r.average.p,
    averages.f = r.average.f,
    baseline.matrix.p = baseline.matrix.p,
    observation.matrix.p = Positives,
    baseline.matrix.f = baseline.matrix.f,
    observation.matrix.f = Fails,
    map = map,
    column = i,
    region_names = LTLA_names,
    n.matrix.p = N.p,
    n.matrix.f = N.f)
  dev.off()
}
