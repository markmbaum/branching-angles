/*
This is a copy of the script used to collect NDVI and SMAP data
in the Google Earth Engine code editor. It provides exports of
 1. maximum NDVI at ~5 km resolution
 2. mean SMAP at ~10 km resoluion
both clipped to the continental United States.
*/

//bounding box for continental Unites States
var conus = ee.Geometry.BBox(-125.48, 24.86, -65.93, 49.84)

//global 16-day, 1 km, NDVI data
var ndvi_max = ee.ImageCollection('MODIS/006/MOD13A2')
  //take maximum values and clip to the CONUS box
  .max().clip(conus)
  //unsure why this combination of steps is necessary,
  //but it reduces the resolution
  .reproject({
    crs: 'EPSG:4326',
    scale: 5000
  })
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    maxPixels: 1024
  })

var smap_mean = ee.ImageCollection("NASA_USDA/HSL/SMAP10KM_soil_moisture").mean().clip(conus)

Map.addLayer(ndvi_max)
Map.addLayer(smap_mean)

Export.image.toDrive({
  image: ndvi_max,
  description: 'conus_ndvi',
  scale: 5000
});

Export.image.toDrive({
  image: smap_mean,
  description: 'smap_mean',
  scale: 10000
});