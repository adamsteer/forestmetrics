#! /usr/bin/python3
"""
Utility functions for generating TERN forest products
from classified lidar point clouds

Dr Adam Steer
November 2019
"""

import pdal
import numpy as np
import json

import shapely
from shapely.geometry import Point
from shapely.geometry import MultiPolygon
from shapely.geometry import box
#from shapely.strtree import STRtree

import geopandas as gpd
import pandas as pd
import osmnx as ox

import os

import sys

# this is needed to create a raster from the output array
from osgeo import gdal
import osgeo.osr as osr
import boto3

def get_elvis_data(elvisresult, aws_params):


    s3 = boto3.client('sts',
                aws_access_key_id=aws_params["aws_key"],
                aws_secret_access_key= aws_params["aws_secret"],
                region_name=aws_params["region"]
                )
    assumed_role_object = s3.assume_role(
                RoleArn=aws_params["aws_role_arn"],
                ExternalId=aws_params["aws_external_id"],
                RoleSessionName=aws_params["aws_role_session_name"]
                )
    credentials = assumed_role_object['Credentials']
    s3_resource = boto3.resource('s3',
                aws_access_key_id=credentials['AccessKeyId'],
                aws_secret_access_key=credentials['SecretAccessKey'],
                aws_session_token=credentials['SessionToken'],
                region_name='ap-southeast-2')

    bucketparts = fileurl.split("/")

    return()

def write_product_geotiff(griddedpoints, outfile, parameters):
    """
    writes out a geotiff from a numpy array of forest metric
    results.

    inputs:
    - a numpy array of metrics [griddedpoints]
    - an outfile name [outfile]
    - a dictionary of parameters for the raster

    outputs:
    - a gdal dataset object
    - [outfile] written to disk
    """

    width = parameters["width"]
    height = parameters["height"]
    wktcrs = parameters["projection"]

    srs = osr.SpatialReference()
    srs.ImportFromWkt(wktcrs)

    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(outfile, width, height, 1, gdal.GDT_Float32 )
    ds.SetGeoTransform([parameters["upperleft_x"],
                       parameters["resolution"],
                       0,
                       parameters["upperleft_y"],
                       0,
                       -parameters["resolution"]])
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(np.rot90(griddedpoints))

    ds.FlushCache()
    ds = None

    return()

#convert PDAL structured array to GeoPandas dataframe
def pdal2df(points):
    """
    Feed me a PDAL pipeline return array, get back a
    GeoPandas dataframe

    Inputs:
    - a PDAL structured array

    Output:
    - a GeoPandas GeoDataFrame with a geometry column generated from
      X and Y coordinates
    """

    arr = points[0]
    #description = arr.dtype.descr
    #cols = [col for col, __ in description]
    gdf = gpd.GeoDataFrame(arr)
    gdf.name = 'nodes'
    #gdf['geometry'] = gdf.apply(lambda row: Point((row['X'], row['Y'])), axis=1)
    #adding a geometry for each point is a slow step...
    gdf["geometry"] = [shapely.geometry.Point(xy) for xy in zip(gdf.X, gdf.Y)]

    gdf = gdf[["X", "Y", "Z", "HeightAboveGround", "Classification", "Intensity", "ReturnNumber", "NumberOfReturns", "geometry"]]

    return(gdf)

# create a GeoPandas rTree spatial index
def spatialindex(dataframe):
    """
    create an rTree index for a GeoPandas GeoDataFrame

    Input:
    - a GeoPandas GeoDataFrame with a geometry column

    Output:
    - a GeoPandas rTree spatial index for the input dataframe
    """
    sindex = dataframe.sindex
    return(sindex)

#get a pointview from PDAL
def readlasfile(lasfile):
    """
    Run a PDAL pipeline. Input is a JSON declaration to
    deliver to PDAL. Output is a labelled numpy array.

    Data are filtered to compute height above ground using nearest ground point neighbours (TIN method arriving soon) and sort by morton order. Any unused dimensions are also trimmed.
    """
    pipeline = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": lasfile
            },
            {
                "type": "filters.hag"
            },
            {
                "type": "filters.mortonorder"
            }
        ]
    }

    #create a pipeline object
    pipeline = pdal.Pipeline(json.dumps(pipeline))

    # execute the pipeline
    count = pipeline.execute()

    #read points into a numpy structured array
    arrays = pipeline.arrays

    #return the numpy array to operate on
    return(arrays)

#read LAS/LAZ metadata
def readlasmetadata(lasfile):
    """
    read only the LAS file metadata for a given tile

    Inputs:
    - a LAS file name

    Outputs:
    - a dictionary object containing LAS file metadata fields

    Notes:
    In this function only the first point in the input LAS/LAZ file is
    read. This provides fast access to LAS metadata, but only gives
    statistical information about the first point.
    To see descriptive stats for a whole LAS tile, remove the 'count'
    parameter from the pipeline JSON.
    Doing so will make metadata reads substantially slower.
    """

    pipeline = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": lasfile,
                "count": 1
            },
            {
                "type": "filters.info"
            }
        ]
    }

    pipeline = pdal.Pipeline(json.dumps(pipeline))
    count = pipeline.execute()

    #extract metadata into a JSON blob
    metadata = json.loads(pipeline.metadata)

    return(metadata)

def gen_raster_cells(metadata, resolution):
    """
    Generate cells of 'resolution x resolution' for point querying

    input:
    - PDAL metadata

    output:
    - shapely geometry containing polygons defining 'resolution x resolution'
      boxes covering the LAS tile extent

    """
    bbox = box(metadata["metadata"]["readers.las"]["minx"],
               metadata["metadata"]["readers.las"]["miny"],
               metadata["metadata"]["readers.las"]["maxx"],
               metadata["metadata"]["readers.las"]["maxy"])

    tiledBBox = ox.quadrat_cut_geometry(bbox, quadrat_width=resolution)

    return(tiledBBox)

def get_cell_points(poly, df, sindex):
    """
    Query a spatial index for points within a shapely polygon

    """
    poly = poly.buffer(1e-14).buffer(0)
    possible_matches_index = list(sindex.intersection(poly.bounds))
    possible_matches = df.iloc[possible_matches_index]
    precise_matches = possible_matches[possible_matches.intersects(poly)]

    return(precise_matches)

def comp_dem(lasfile, outpath, resolution):
    # interpolate ground returns in a grid and output a raster
    # this is likely to be handled by PDAL... or GDAL

    fileroot = make_file_rootname(lasfile)

    outfilename = os.path.join(outpath,
                               "dem",
                               fileroot + "-DEM-" + str(resolution) + "m.tif")

    pipeline = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": lasfile
            },
            {
                "type": "filters.range",
                "limits": "Classification[2:2]",
            },
            {
                "type": "writers.gdal",
                "filename": outfilename,
                "resolution": resolution,
                "output_type": "idw",
                "window_size": 3
            }
        ]
    }

    #create a pipeline object
    pipeline = pdal.Pipeline(json.dumps(pipeline))
    count = pipeline.execute()

    return()

def make_file_rootname(lasfile):
    """
    make a file 'root name' from an input file path

    input:
    - a path to an input file (string)

    output:
    - a string containing the file name minus extension, assuming
    the extension is 4 characters long
    """
    filebits = lasfile.split("/")
    infilename = filebits[-1]
    fileroot = infilename[:-4]
    return(fileroot)

def read_data(lasfile):
    """
    #wrapper to read in LAS data and produce a dataframe + spatial index
    """
    metadata = readlasmetadata(lasfile)

    points = readlasfile(lasfile)

    dataframe = pdal2df(points)

    spatial_index = spatialindex(dataframe)

    return(metadata, dataframe, spatial_index)

def compute_tern_products(metadata, points, sindex, resolution, lasfile, outpath):
    """
    Wrapper to iterate over the input data and generate rasters for each product.

    *note this part could be paralellised - maybe per-product, or per-cell

    Each grid square processed in this loop corresponds to one pixel in an output raster.

    """

    #set up an 'output resolution' sized grid - like a fishnet grid.
    # each polygon in the resulting set covers an area of 'resolution X resolution'
    pixel_grid = forestutils.gen_raster_cells(metadata, resolution)

    #set up output rasters

    # get tile width and height
    tile_width = metadata["metadata"]["readers.las"]["maxx"] - metadata["metadata"]["readers.las"]["minx"]
    tile_height = metadata["metadata"]["readers.las"]["maxy"] - metadata["metadata"]["readers.las"]["miny"]

    raster_xsize = int(np.ceil(tile_width) / resolution)
    raster_ysize = int(np.ceil(tile_height) / resolution)

    #replicate for all products...
    vh_raster = np.zeros((raster_xsize, raster_ysize))
    vcf_raster = np.zeros((raster_xsize, raster_ysize))
    cth_raster = np.zeros((raster_xsize, raster_ysize))
    cbh_raster = np.zeros((raster_xsize, raster_ysize))
    fbf_raster = np.zeros((raster_xsize, raster_ysize))
    cli_raster = np.zeros((raster_xsize, raster_ysize))
    density_raster = np.zeros((raster_xsize, raster_ysize))

    veg_below_dict = {}

    veg_below_dict["all"] = np.zeros((raster_xsize, raster_ysize))
    for height in LCF_HEIGHTS:
        veg_below_dict[str(height)] = np.zeros((raster_xsize, raster_ysize))

    #internal loop around grid squares covering the LAS tile.
    # this is another ppoint for parallelisation - since we can set up a list of geometries
    # and cast that at multipuple processes, setting up one process per grid square
    # another way to do this would be to recast this loop block into a function which can
    # be called by one process per product
    # the second strategy seems easier, then only one process is trying to write into each
    # output array.

    for pixel in pixel_grid:

        #compute output array index for this cell:
        poly_x, poly_y = pixel.centroid.xy

        poly_base_x = poly_x[0] - metadata["metadata"]["readers.las"]["minx"]
        poly_base_y = poly_y[0] - metadata["metadata"]["readers.las"]["miny"]

        array_x = int(np.floor((poly_base_x / (resolution)) ))
        array_y = int(np.floor((poly_base_y / (resolution)) ))

        #get points for this cell
        matches = forestutils.get_cell_points(pixel, points, sindex)

        #compute in order
        #VH

        vh_raster[array_x, array_y] = metrics.comp_vh(matches)

        #VCF
        vcf_raster[array_x, array_y] = metrics.comp_vcf(matches)

        #LCF - long-ish process..
        # compute a dictionary of points below height thresholds
        veg_below = metrics.comp_veg_layers(matches, LCF_HEIGHTS)

        # add the first element of the dictionary to a raster output
        veg_below_dict["all"][array_x, array_y] = veg_below["all"]

        #iterate over the height thresholds and do likewise...
        for height in LCF_HEIGHTS:
            veg_below_dict[str(height)][array_x, array_y] = veg_below[str(height)]

        #CTH
        cth_raster[array_x, array_y] = metrics.comp_cth(matches)

        #CBH
        cbh_raster[array_x, array_y] = metrics.comp_cbh(matches)

        #FBF
        fbf_raster[array_x, array_y] = metrics.comp_fbf(matches)

        #CLI
        cli_raster[array_x, array_y] = metrics.comp_cli(matches)

        #density
        density_raster[array_x, array_y] = metrics.comp_density(matches, resolution)

    #end of computing stuff, time to make outputs...

    #compute LCF values given our height thresholded veg counts
    lcf = metrics.comp_lcf(veg_below_dict, vcf_raster)

    if (not os.path.isdir(outpath + "/dem")):
        os.mkdir(outpath + "/dem")

    dem = forestutils.comp_dem(lasfile, outpath, resolution)

    tern_products = {}
    tern_products["vh"] = vh_raster
    tern_products["vcf"] = vcf_raster
    tern_products["cth"] = cth_raster
    tern_products["cbh"] = cbh_raster
    tern_products["fbf"] = fbf_raster
    tern_products["cli"] = cli_raster
    tern_products["lcf_h"] = lcf["lcf_h"]
    tern_products["lcf_os"] = lcf["lcf_os"]
    tern_products["lcf_us"] = lcf["lcf_us"]
    tern_products["lcf_cf"] = lcf["lcf_cf"]
    tern_products["lcf_ef"] = lcf["lcf_ef"]
    tern_products["lcf_nsf"] = lcf["lcf_nsf"]
    tern_products["density"] = density_raster

    return(tern_products)

def export_tern_products(tern_products, metadata, resolution, lasfile, outpath):

    #set up GDAL parameters

    wktcrs = metadata["metadata"]["readers.las"]["comp_spatialreference"]

    raster_parameters = {}
    raster_parameters["width"] = np.shape(tern_products["vh"])[0]
    raster_parameters["height"] = np.shape(tern_products["vh"])[1]
    raster_parameters["upperleft_x"] = metadata["metadata"]["readers.las"]["minx"]
    raster_parameters["upperleft_y"] = metadata["metadata"]["readers.las"]["maxy"]
    raster_parameters["resolution"] = resolution
    raster_parameters["projection"] = wktcrs

    fileroot = forestutils.make_file_rootname(lasfile)
    print(fileroot)

    for productname in tern_products.keys():

        if (not os.path.isdir(os.path.join(outpath,
                                       productname))):
            os.mkdir(os.path.join(outpath,
                                  productname))


        #set output filenames
        separator = "-"

        raster_name = separator.join([fileroot,
                                     productname,
                                     str(resolution) + "m.tiff"])
        raster_path = os.path.join(outpath,
                                   productname,
                                   raster_name)
        print(raster_path)
        forestutils.write_product_geotiff(tern_products[productname], raster_path, raster_parameters)


    return()
