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
    description = arr.dtype.descr
    cols = [col for col, __ in description]
    gdf = gpd.GeoDataFrame({col: arr[col] for col in cols})
    gdf.name = 'nodes'
    gdf['geometry'] = gdf.apply(lambda row: Point((row['X'], row['Y'], row['Z'])), axis=1)

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

    Data are filtered to compute height above ground using nearest ground point neighbours (TIN method arriving soon)
    """
    pipeline = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": lasfile
            },
            {
                "type": "filters.hag"
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

    # execute the pipeline
    try:
        count = pipeline.execute()
    except ValueError:
        print("DEM could not be written")

    return()

def make_file_rootname(lasfile):
    filebits = lasfile.split("/")
    infilename = filebits[-1]
    fileroot = infilename[:-4]
    return(fileroot)

def read_data(lasfile):
    """
    wrapper to read in LAS data and produce a dataframe + spatial index
    """

    metadata = readlasmetadata(lasfile)

    points = readlasfile(lasfile)

    dataframe = pdal2df(points)

    spatial_index = spatialindex(dataframe)

    return(metadata, dataframe, spatial_index)
