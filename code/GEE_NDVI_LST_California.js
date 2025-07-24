title: "Quantifying the Impacts of Building Energy Efficiency Retrofits and Nature Exposure on Chronic Student Absenteeism"
author: "Devan Cantrell Addison-Turner, Yingjie Li, Anthony Dylan Kinslow II, Gretchen Cara Daily, and Rishee Kumar Jain"

// Code to calculate the Normalized Difference Vegetation Index (NDVI) and Land Surface Temperature (LST) near K-12 Schools 
// across California, USA from satellite imagery using Google Earth Engine (GEE) during the COVID-19 Pandemic.

var dir_output = "Devan_gee_data"
var version = 'v3_2'

// Define the time range
var startYear = 2020;
var endYear = 2022;

// shapefile for CA
var tracts = ee.FeatureCollection('TIGER/2020/TRACT')
      .filter(ee.Filter.eq("STATEFP", "06"))


var tracts = ee.FeatureCollection("TIGER/2018/States")
      .filter(ee.Filter.eq("STATEFP", "06"))
    
// var school = ee.FeatureCollection("users/daddisonturner/CA_Schools_600_meter_buffer_GIS") // used before 6/17
    // .select(['id', 'CDS_Code', 'School_ID'])
var school = ee.FeatureCollection("users/daddisonturner/SCE_impacts_V3_shapefiles")
    .select(['id', 'CDSCode', 'UID'])
    // .limit(100) // subset the first 100 rows for testing
var geoid = "School_ID"
var geoid = "UID"

// var geometry = tracts
// var geometry = school_filtered; // for testing use
var geometry = school; 
// print('geometry', geometry);
Map.addLayer(geometry, {}, 'geometry');


// ##########################################################################################################
// Download NDVI image using Landsat 8 for any study region 
// ##########################################################################################################

// Cloud masking function for Landsat 8 Collection 2 Level 2
function maskL8sr(image) {
  var cloudShadowBitMask = (1 << 4);
  var cloudsBitMask = (1 << 3);
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                   .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask)
              .copyProperties(image, ['system:time_start']);
}

var applyScaleFactors = function(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBand = image.select('ST_B10').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBand, null, true);
}


// 1. Import the Landsat 8 image 
var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
// 2. Get the least cloudy image in 2020 - 2022
    // .filterBounds(geometry)
    .filterDate('2020-01-01', '2023-01-01')
    .map(applyScaleFactors)
    .map(maskL8sr)
    .map(function(image){return image.clip(geometry)});

// 3. Compute the Normalized Difference Vegetation Index (NDVI)
var L8_ndvi = L8
      .median()
      .normalizedDifference(['SR_B5', 'SR_B4'])
      .rename('NDVI')
      .select('NDVI')
// print('L8_ndvi', L8_ndvi)

// 4. Display the result
var palette = ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718',
               '74A901', '66A000', '529400', '3E8601', '207401', '056201',
               '004C00', '023B01', '012E01', '011D01', '011301'];
// print('L8_ndvi', L8_ndvi)
Map.addLayer(L8_ndvi, {min: 0, max: 1, palette: palette},'NDVI'); 

//// 5. Export to Drive
// Export.image.toDrive({
//   image: L8_ndvi,
//   description: 'NDVIimage_2020_2022',
//   folder: 'California NDVI',
//   scale: 30,
//   region: geometry,
//   maxPixels: 1e13,
// });




// ##########################################################################################################
// Statistics of Image Regions
// ##########################################################################################################

// You can combine reducers to calculate e.g. mean and standard deviation
// simultaneously. The resulting property names are the concatenation of the
// band names and statistic names, separated by an underscore.
var reducer = ee.Reducer.mean().combine({
  reducer2: ee.Reducer.median(),
  sharedInputs: true
});


// Add reducer output to the Features in the collection.
var ndvi_mean = L8_ndvi.reduceRegions({
  collection: geometry,
  reducer: reducer,
  scale: 30,
});


// Create an empty image into which to paint the features, cast to byte.
var empty = ee.Image().byte();
// Paint the edges with different colors, display.
var ndvi_viz = empty.paint({
  featureCollection: ndvi_mean,
  color: 'mean',
  // width: 1 // If the width parameter is not provided, the interior of the features is painted
})
Map.addLayer(ndvi_viz, {palette: palette, max: 1}, 'ndvi_viz mean', false);

var ndvi_viz = empty.paint({
  featureCollection: ndvi_mean,
  color: 'median'
})
Map.addLayer(ndvi_viz, {palette: palette, max: 1}, 'ndvi_viz median', false);


// export results

// for l8
var DataRenamed = ndvi_mean.map(function(feat){
  return ee.Feature(feat.geometry(), { 
    ndvi_mean: feat.get('mean'),
    ndvi_median: feat.get('median'),
    GEOID: feat.get(geoid)
  })
})

var filename = 'ndvi_by_regions_' + 'CA_2020_' + geoid + '_' + version; print(filename)
Export.table.toDrive({
           collection: DataRenamed, 
           description: filename, 
           fileFormat: 'CSV', 
           folder: dir_output, // result folder
           selectors: ['GEOID', 'ndvi_mean', 'ndvi_median']})
           








/// ######################################################################################################################         
/// Land Surface Temperature 
/// ######################################################################################################################  
// ref: https://gis.stackexchange.com/questions/422033/land-surface-temperature-in-google-earth-engine
// ref: https://medium.com/@ridhomuh002/analyzing-land-surface-temperature-lst-with-landsat-8-data-in-google-earth-engine-f4dd7ca28e70


//vis params
var vizParams2 = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0,
  max: 3000,
  gamma: 1.4,
};

//load the collection:
var col = L8
//imagen reduction
//median
// var image = col.median();
// print('image', image);
// Map.addLayer(image, vizParams2, 'image');

// Function to filter collection by year and calculate median
function getYearlyMedian(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year, 12, 31);
  
  var yearCollection = col.filterDate(startDate, endDate)
                                        .filterBounds(geometry)
                                        .map(maskL8sr)
                                        .median();
                                        
  return yearCollection.set('year', year);
}

// Create a list of years
var years = ee.List.sequence(startYear, endYear);

// Calculate the median for each year
var col_yearlyMedians = ee.ImageCollection(years.map(getYearlyMedian));
// // Get the size of the original Landsat collection for the specified date range and region
// var collectionSize = col_yearlyMedians.size();
// print('Size of the Collection:', collectionSize);


//individual LST images
var col_list = col_yearlyMedians.toList(col.size());
// print('size', col_yearlyMedians.size());


var LST_col = col_list.map(function (ele) {
  var date = ee.Image(ele).get('system:time_start');
  var ndvi = ee.Image(ele).normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  // find the min and max of NDVI
  var min = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: geometry,
    scale: 30,
    maxPixels: 1e9
  }).values().get(0));
  
  var max = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: geometry,
    scale: 30,
    maxPixels: 1e9
  }).values().get(0));
  
  
  // Fraction of Vegetation (FV) Calculation
  // Formula: ((NDVI - NDVI_min) / (NDVI_max - NDVI_min))^2
  // Calculate the Fraction of Vegetation (FV) using the NDVI values within the specified range.
  // NDVI_min represents the minimum NDVI value, and NDVI_max represents the maximum NDVI value
  var fv = (ndvi.subtract(min).divide(max.subtract(min)))
              .pow(ee.Number(2))
              .rename('FV');
  
  
  // Emissivity Calculation
  // Formula: 0.004 * FV + 0.986
  // Calculate Land Surface Emissivity (EM) using the Fraction of Vegetation (FV).
  // The 0.004 coefficient represents the emissivity variation due to vegetation,
  // and the 0.986 represents the base emissivity for other surfaces.
  var a= ee.Number(0.004);
  var b= ee.Number(0.986);
  
  var EM = fv.multiply(a).add(b).rename('EM').select('EM');

  var image = ee.Image(ele);
  // Select Thermal Band (Band 10) and Rename It
  var thermal = image.select('ST_B10').rename('thermal');
  
  
  // Now, lets calculate the land surface temperature (LST)
  // Formula: (TB / (1 + (λ * (TB / 1.438)) * ln(em))) - 273.15
  var LST = image.expression(
    '(TB/(1 + (0.00115* (TB / 1.438))*log(em)))-273.15', {
      'TB': thermal.select('thermal'), // Select the thermal band (TB),
      'em': EM // Assign emissivity (em)
  });

  return ee.Algorithms.If(min, LST.set('system:time_start', date).float().rename('LST'), 0);
  // return LST.set('system:time_start', date);
}).removeAll([0]);

LST_col = ee.ImageCollection(LST_col);



////--------------------------------------
//// --- one image as input --------------
////--------------------------------------
var ele = L8.median().clip(geometry);
var ndvi = ee.Image(ele).normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
// find the min and max of NDVI
var min = ee.Number(ndvi.reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: geometry,
  scale: 30,
  maxPixels: 1e9
}).values().get(0));

var max = ee.Number(ndvi.reduceRegion({
  reducer: ee.Reducer.max(),
  geometry: geometry,
  scale: 30,
  maxPixels: 1e9
}).values().get(0));
  
var fv = (ndvi.subtract(min).divide(max.subtract(min)))
            .pow(ee.Number(2))
            .rename('FV');

var a= ee.Number(0.004);
var b= ee.Number(0.986);
var EM = fv.multiply(a).add(b).rename('EM');
var image = ee.Image(ele);
var thermal = image.select('ST_B10').rename('thermal');
var LST = image.expression(
  '(TB/(1 + (0.00115* (TB / 1.438))*log(em)))-273.15', {
    'TB': thermal.select('thermal'), // Select the thermal band (TB),
    'em': EM // Assign emissivity (em)
}).float().rename('LST');

Map.addLayer(LST, {min: 0, max: 50, palette: [
    '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
    '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
    '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
    'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
    'ff0000', 'de0101', 'c21301', 'a71001', '911003']},'LST');





//print(LST);
// var LST_Params = {min: 25, max: 40, palette: ['blue', 'white', 'red']};
var LST_Params = {
  min: 0, // Minimum LST value
  max: 50, // Maximum LST value
  palette: [
    '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
    '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
    '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
    'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
    'ff0000', 'de0101', 'c21301', 'a71001', '911003'
  ]}

// var LST_mean = LST_col.mean();
// Map.addLayer(LST_mean, LST_Params, 'LST_mean');

// var export_Image = LST_col.reduce(ee.Reducer.mean());
var export_Image = LST.clip(geometry)
// Map.addLayer(export_Image, LST_Params, 'LST_mean export_Image');
// Map.addLayer(LST, LST_Params, 'LST');





var crs_Projection = export_Image.projection();
// Get the forest cover data at MODIS scale and projection.
var export_Image_resample = export_Image
    // Force the next reprojection to aggregate instead of resampling.
    // .reduceResolution({
    //   reducer: ee.Reducer.mean(),
    //   maxPixels: 1024
    // })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: crs_Projection,
      scale: 1000
    });
// print("export_Image", export_Image);
// print("export_Image_resample", export_Image_resample);


// To summarize the image by school regions and export the results as a csv file
// ######################################
// Statistics of Image Regions
// ######################################

// Add reducer output to the Features in the collection.
var mean_byRegions = export_Image.reduceRegions({
  collection: geometry,
  reducer: reducer,
  scale: 30,
});
// print('mean_byRegions', mean_byRegions) 
/// ! do not print if using geometry of all schools -> take forever and result in errors


// // Create an empty image into which to paint the features, cast to byte.
// var empty = ee.Image().byte();
// // Paint the edges with different colors, display.
// var res_viz = empty.paint({
//   featureCollection: mean_byRegions,
//   color: 'mean',
//   // width: 1 // If the width parameter is not provided, the interior of the features is painted
// })
// // Map.addLayer(res_viz, LST_Params, 'res_viz mean', false);

// var res_viz = empty.paint({
//   featureCollection: mean_byRegions,
//   color: 'median'
// })
// Map.addLayer(res_viz, LST_Params, 'res_viz median', false);




// ######################################
// export results
// ######################################

var DataRenamed = mean_byRegions.map(function(feat){
  return ee.Feature(feat.geometry(), { 
    mean_byRegions:   feat.get('mean'),
    median_byRegions: feat.get('median'),
    GEOID: feat.get(geoid)
  })
})

var filename = 'LST_byRegions_' + 'CA_2020_' + geoid + '_' + version; print(filename)
Export.table.toDrive({
           collection: DataRenamed, 
           description: filename, 
           fileFormat: 'CSV', 
           folder: dir_output, // result folder
           selectors: ['GEOID', 'mean_byRegions', 'median_byRegions']})
