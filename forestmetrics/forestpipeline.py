#! /usr/bin/python3
"""
Pipeline for running forest structure metrics



"""

import forestmetrics
import utils

import numpy as np

import pdal
from osgeo import gdal
import shapely

import boto3
import os
import sys
import zipfile
# callingelvise needs to be on the Python path
#from callingelvis.anybodyhome import anybodyhome

### using a list of file names...
def forestpipeline(input_lidar):

    for lidar_file in input_lidar:
        if (lidar_file[:-4].find(zip)):
            ###unzip the file...
            zipfile.Zipfile.extract(lidar_file)
            #replace .zip with .las in variable name
            lidarfile.replace(".zip", ".las")

        thepoints = utils.readlasfile(lidar_file)
        themetadata = utils.readlasmetadata(lidar_file)

        pointsgdf = utils.pdal2df(thepoints)
        pointsindex = utils.spatialindex(pointsgdf)

"""
elvisparams = []

bbox = elvisparams["bbox"]
pointtype = elvisparams["pointtype"]
collection = elvisparams["collection"]
year = elvisparams["year"]
heightref = elvisparams["heightref"]
"""
## first call elvis
datasets = callelvis(bbox, pointtype, collection, year, heightref)


## next parse results into an array that can be delivered to an
# iterative *or* parallel processing pipeline
lidar_files = datasets[jurisdiction][height]["filename"]

## next fetch results and process.
