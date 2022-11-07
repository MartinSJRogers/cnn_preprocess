# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 14:02:44 2022

@author: marrog

Rasterise shapefile using corresponding rasrer image
Raster image could be visible or SAR image. Use this image to get important
information e.g. resolution and projection so that it is possible to export
the rasterised shapefile as a geotiff. 
"""

from osgeo import gdal
from osgeo import ogr
from osgeo import gdalconst
import geopandas
import rasterio


# rasterise polygon shapefiles at 40m resolution. Export to same directory
# as the reference images dataset
def rasteriseShapefile():

    # manually import file names for SAR image and polygonised shapefile
    sar_fn = "*.tif"
    shp_fn = "*.shp"
    # define output rasterised shapefile name
    output = "*.tif"

    # open shapefile and corresponding raster
    data = gdal.Open(sar_fn, gdalconst.GA_ReadOnly)
    dataArray = data.ReadAsArray()
    raster = rasterio.open(sar_fn)

    # get raster epsg to convert polygon shapefile into
    crsOut = raster.crs
    epsgOut = str(crsOut)[-4:]
    print(epsgOut)

    # Get coordinate information for input raster
    # x_min/x_max, y_min/y_max correspond to bounding box of image
    # x and y res are the verticle and horizontal spatial resolution
    geo_transform = data.GetGeoTransform()
    x_min = geo_transform[0]
    y_max = geo_transform[3]
    x_max = x_min + geo_transform[1] * data.RasterXSize
    y_min = y_max + geo_transform[5] * data.RasterYSize
    x_res = data.RasterXSize
    y_res = data.RasterYSize
    pixel_width = geo_transform[1]

    # Open shapefile in geopandas
    dataShp = geopandas.read_file(shp_fn)
    # change CRS from polygon to raster file
    dataShp = dataShp.to_crs(epsg=epsgOut)
    # convert geopandas into ogr Datasource
    shp_temp = ogr.Open(dataShp.to_json())
    # extract layer from ogr Datasource
    mb_l = shp_temp.GetLayer()

    # create dataset as geotiff
    target_ds = gdal.GetDriverByName('GTiff').Create(
        output, x_res, y_res, 1, gdal.GDT_Byte)
    # attach geo information to dataset, e.g. resolution and projection
    target_ds.SetGeoTransform((x_min, pixel_width, 0, y_min, 0, pixel_width))
    band = target_ds.GetRasterBand(1)
    NoData_value = -999999
    band.SetNoDataValue(NoData_value)
    band.FlushCache()
    # Need to change the value after 'ATTRIBUTE=". This corresponds to the field (column)
    # in the shapefile's attribute table that you want to use.
    gdal.RasterizeLayer(target_ds, [1], mb_l, options=["ATTRIBUTE=xxx"])

    # important to do this to precent errors coming up
    target_ds = None


rasteriseShapefile()
