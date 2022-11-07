# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 15:41:59 2021

@author: marrog

Main preprocessing in CNN pipeline
Place all SAR images ande corresponding shapefiles into list alongside data

rasterise all shapefiles using same resolution as original SAR image

Use sliding window to crop the large images and corresponding rasterised shapefiles
    to 720 by 720 pixels (or size of your choice). Place all images into list.
"""
from osgeo import gdal
from osgeo import ogr
from osgeo import gdalconst
import numpy as np
import rasterio
import geopandas
import os
import glob
from PIL import Image


class Preprocess():
    def __init__(self):
        # manually predefined parameters, can edit and set up using argParse etc.
        self.patchHeight = 720
        self.patchWidth = 720
        self.overlap = 60
        # suggested directory structure
        self.root_dir = "*/*/images"
        self.postprocess_dir = "*/*/images/postprocess"
        self.output_dir = "*/*/images/postprocess/patches"

    # Extract all SAR images and manually digitised boundaries (polygon90) from subdirecties
    # the if == True command ensures images only from directories with both images is collected
    # this is specific to me, but I ensured that original SAR images and corresponding shapefiles
    # ended with the same string, so that I could place them in a list.
    def getFileNames(self):

        path = self.root_dir
        resultRaster = []
        resultPoly = []
        dateTime = []
        imgNumber = []

        os.chdir(path)
        i = 100

        for fileRaster in glob.glob("**/ASAR_WSM.dim.tif", recursive=True):
            subdirectory = os.path.normpath(
                str(self.root_dir) + '/' + str(fileRaster[:-17]))
            for File in os.listdir(subdirectory):
                if File.endswith("polygon90.shp"):

                    i = i+1
                    imgNumber.append(i)
                    # make list of all dates (so these can be used in filenames for output images)
                    dateTime.append(fileRaster[15:30])
                    dir_name = os.path.dirname(os.path.abspath(fileRaster))
                    # make list of all polygons
                    resultPoly.append(str(dir_name) + '/polygon90.shp')
                    # make list of all raster images
                    resultRaster.append(os.path.join(path, fileRaster))

        return resultRaster, resultPoly, dateTime, imgNumber


# rasterise polygon shapefiles at resolution of original SAR image. Export to same directory
# as the reference images dataset


    def rasterise_multi(self, polyList, rasterList, dateTime, imgNumber):
        path = self.root_dir
        outputDirectory = self.postprocess_dir
        date = dateTime
        # concatenate lists so you have (polygon file name, raster file name, date of image, file number (arbitary))
        filesComb = list(zip(polyList, rasterList, dateTime, imgNumber))

        # Run through concatenated list. Rasterise each polygon layer using same code as in
        # rasterise_shapefiles.py code.
        for i in filesComb:
            # can use this if statement so that only every 100th image in list ( or other number),
            # is analysed to check the process is runnign correctly
            # if i[3]%100==0:
            shp1 = str(i[0])
            ndsm = str(i[1])
            date = str(i[2])

            # open shapefile and raster

            data = gdal.Open(ndsm, gdalconst.GA_ReadOnly)
            dataArray = data.ReadAsArray()
            raster = rasterio.open(ndsm)
            # get raster epsg to convert polygon shapefile into
            crsOut = raster.crs
            epsgOut = str(crsOut)[-4:]

            # define output shapefile name
            output = str(outputDirectory) + '/' + str(date) + 'ref.tif'

            # Get coordinate information for input raster
            geo_transform = data.GetGeoTransform()

            #source_layer = data.GetLayer()
            x_min = geo_transform[0]
            y_max = geo_transform[3]
            x_max = x_min + geo_transform[1] * data.RasterXSize
            y_min = y_max + geo_transform[5] * data.RasterYSize  # (-40)
            x_res = data.RasterXSize
            y_res = data.RasterYSize
            pixel_width = geo_transform[1]

            # Open shapefile in geopandas
            dataShp = geopandas.read_file(shp1)
            # change CRS from polygon to raster file
            dataShp = dataShp.to_crs(epsg=3031)
            # convert geopandas into ogr Datasource
            shp_temp = ogr.Open(dataShp.to_json())
            # extract layer from ogr Datasource
            mb_l = shp_temp.GetLayer()

            target_ds = gdal.GetDriverByName('GTiff').Create(
                output, x_res, y_res, 1, gdal.GDT_Byte)

            target_ds.SetGeoTransform(
                (x_min, pixel_width, 0, y_min, 0, pixel_width))
            band = target_ds.GetRasterBand(1)
            NoData_value = -999999
            band.SetNoDataValue(NoData_value)
            band.FlushCache()
            # again change attribute name if needed
            gdal.RasterizeLayer(
                target_ds, [1], mb_l, options=["ATTRIBUTE=type"])

            target_ds = None

            # normalise input SAR image to [-1 , 1] prior to saving
            normalized_input = (dataArray - np.amin(dataArray)) / \
                (np.amax(dataArray) - np.amin(dataArray))
            normalised = (2*normalized_input) - 1

            im = Image.fromarray(normalised)
            im.save(str(outputDirectory) + '/' + str(date) + 'raw.tif')

        return filesComb

    # Use sliding window to crop the large images and corresponding rasterised shapefiles
    # to 720 by 720 pixels (or size of your choice). Place all images into list.
    def save_pan_patches(self):

        # extract user defined information from _init_ function
        preprocess_path = self.postprocess_dir
        save_path = self.output_dir
        patch_height = self.patchHeight
        patch_width = self.patchWidth
        overlap = self.overlap

        # create empty numpy arrays for the saved image names to be listed ion
        allNames = np.empty([1, 2])
        edgeNames = np.empty([1, 2])
        iceNames = np.empty([1, 2])
        waterNames = np.empty([1, 2])
        equalNames = np.empty([1, 2])

        os.chdir(preprocess_path)
        # run through every SAR image...
        for fileRaster in glob.glob("**/*raw.tif", recursive=True):
            # extract file name for original SAR image (rawImg) and rasterised polygon (refImg)
            rawImg = os.path.normpath(
                str(preprocess_path) + '/' + str(fileRaster))
            refImg = os.path.normpath(
                str(preprocess_path) + '/' + str(fileRaster[:-7]) + 'ref.tif')

            # extract string of date of image capture (used for output file name)
            date = str(fileRaster[-22:-7])

            # convert images into gdal and numpy objects
            imageRef = gdal.Open(refImg)
            arrayRef = imageRef.ReadAsArray()
            arrayRef = np.flip(arrayRef, axis=0)
            imageRaw = gdal.Open(rawImg)
            arrayRaw = imageRaw.ReadAsArray()

            # collect reference image dimensions
            imheight, imwidth = arrayRef.shape

            # n_channels, imheight, imwidth = arrayRaw.shape # when two band SAR (HH and HV)
            imheight, imwidth = arrayRaw.shape

            # set size of stride for cropping kernel
            dx = int(patch_width - overlap*2)
            dy = int(patch_height - overlap*2)
            #print('Slice interval dx = {}, dy = {}'.format(dx, dy))

            n_ims = 0
            # loop through image and crop by 720*720 pixels with a border (overlap) of 120
            # The two if statements ensure that if the 720 pixel crop overlaps the right or
            # bottom of the image, it moves up or left so that it crops the edge of the image
            # and there are no nan values in the crop
            for y in range(0, imheight, dy):
                for x in range(0, imwidth, dx):
                    n_ims += 1

                    if y + patch_height > imheight:
                        y0 = imheight - patch_height
                    else:
                        y0 = y
                    if x + patch_width > imwidth:
                        x0 = imwidth - patch_width
                    else:
                        x0 = x

                    # create numpy array of crop of original SAR (pan_patchRaw) and corresponding
                    # rasterised polygon (pan_patchRef)

                    # pan_patchRaw = arrayRaw[:, y0:y0 + patch_height, x0:x0 + patch_width] # when two band SAR (HH and HV)
                    pan_patchRaw = arrayRaw[y0:y0 +
                                            patch_height, x0:x0 + patch_width]
                    pan_patchRef = arrayRef[y0:y0 +
                                            patch_height, x0:x0 + patch_width]

                    # Only save images which contain no no_data pixels, i.e. that overlap edge of image
                    minVal = np.amin(pan_patchRef)
                    if minVal == 0:
                        print('{}/Patch{}_{}_{}_{}'.format(save_path,
                              n_ims, y0, x0, 'contains no data'))
                    else:
                        # This may not be necessary if you give polygon layers the appropriate values from the start.
                        # Otherwise this npo.where() code can be used to change the pixel values of the rasterised shapefile.

                        # Change values. Currently 0=ND, 1=water 2=ice 9=land 10=unclassified
                        # To: 0=water, 1=ice, 2=land, 3=unclassified, 9999=ND
                        rs1 = np.where(pan_patchRef == 0, 9999, pan_patchRef)
                        rs2 = np.where(rs1 == 1, 0, rs1)
                        rs3 = np.where(rs2 == 2, 1, rs2)
                        rs4 = np.where(rs3 == 9, 2, rs3)
                        pan_patchRef_reVal = np.where(rs4 == 10, 3, rs4)

                        # save image pairs as .npy, can also save as .tif
                        np.save('{}/Patch{}_{}_{}_{}_{}'.format(save_path,
                                date, n_ims, y0, x0, '_Ref'), pan_patchRef_reVal)
                        np.save('{}/Patch{}_{}_{}_{}_{}'.format(save_path,
                                date, n_ims, y0, x0, '_Raw'), pan_patchRaw)

                        refName = 'Patch{}_{}_{}_{}_{}'.format(
                            date, n_ims, y0, x0, '_Ref.npy')
                        rawName = 'Patch{}_{}_{}_{}_{}'.format(
                            date, n_ims, y0, x0, '_Raw.npy')

                        # stack all image names into numpy array.
                        # allNames = every patch
                        # waterNames = patch contains water (at least one pixel = 0)
                        # iceNames = patch contains ice (at least one pixel = 1)
                        newRow = np.array([rawName, refName])
                        allNames = np.vstack([allNames, newRow])
                        if len(np.unique(pan_patchRef_reVal)) > 1:
                            edgeNames = np.vstack([edgeNames, newRow])
                        elif np.max(pan_patchRef_reVal) == 0:
                            waterNames = np.vstack([waterNames, newRow])
                        elif np.max(pan_patchRef_reVal) == 1:
                            iceNames = np.vstack([iceNames, newRow])

        print('Number of patches = {}, sliceheight = {}, slicewidth = {}'.format(
            n_ims, patch_height, patch_width))
        # save file containing every patch name (allNames) and only patches containing both ice and water (edgeNames)
        np.savetxt(r'*/*/allNames.txt',
                   allNames, delimiter='\t', fmt='%s')
        np.savetxt(r'*/*/edgeNames.txt',
                   edgeNames, delimiter='\t', fmt='%s')

        edgeCount = edgeNames.shape[0]
        waterCount = waterNames.shape[0]
        iceCount = iceNames.shape[0]

        shuffleWater = np.random.permutation(waterNames[1:, :])
        shuffleIce = np.random.permutation(iceNames[1:, :])

        # count number of patches that are all water or all ice (every pixel = 0 or 1)
        if waterCount > edgeCount:
            waterEdge = np.vstack(
                [edgeNames[1:, :], shuffleWater[:edgeCount, :]])
        else:
            waterEdge = np.vstack([edgeNames[1:, :], shuffleWater])

        if iceCount > edgeCount:
            equalNames = np.vstack([waterEdge, shuffleIce[:edgeCount, :]])
        else:
            equalNames = np.vstack([waterEdge, shuffleIce])

        # create edgeNames- equal number of all water, all ice and edge patches
        np.savetxt(r'*/*/equalNames.txt',
                   equalNames, delimiter='\t', fmt='%s')

        return iceNames, waterNames, shuffleIce, shuffleWater


a = Preprocess()
resultRaster, resultPoly, dateTime, imgNumber = a.getFileNames()
fileNames = a.rasterise_multi(resultPoly, resultRaster, dateTime, imgNumber)
iceNames, waterNames, shuffleIce, shuffleWater = a.save_pan_patches()
