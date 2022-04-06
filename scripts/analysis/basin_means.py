from pandas import read_csv
from geopandas import read_file

from setup import *

# setup CONUS data
df = read_csv(fnconus)
rename_cols(df)
add_derived_cols(df)

# basin groups and corresponding shapefiles
basins = ['huc2', 'huc4', 'huc6', 'huc8']
fnshps = ['wbdhu2.shp', 'wbdhu4.shp', 'wbdhu6.shp', 'wbdhu8.shp']
# output directory
dirout = join(datadir, 'exp_pro', 'basin-means')

for (basin, fnshp) in zip(basins, fnshps):
    print(basin)

    # average all variables in watersheds and write to csv
    df[basin] = df[basin].astype(int)
    g = df.groupby(basin).mean()
    fn = join(dirout, basin + '_means.csv')
    g.to_csv(fn)
    print('  file written:', fn)

    # join the averages with basins in shapefile and write joined file
    shp = read_file(join(datadir, 'exp_raw', 'us-watersheds', fnshp))
    shp[basin] = shp[basin].astype(int)
    shp = shp.merge(g, on=basin)
    fn = join(dirout, basin + '_means.shp')
    shp.to_file(fn)
    print('  file written:', fn)
