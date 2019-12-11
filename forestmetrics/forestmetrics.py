#! /usr/bin/python3
"""
TERN forest structure metrics

Dr Adam Steer
November 2019

Using the TERN airborne products manual, and MATLAB code produced
by Prof. Albert VanDijk as sources.

"""

#imports

import numpy as np
import json

from shapely.geometry import Point
from shapely.geometry import MultiPolygon
from shapely.geometry import box

#PDAL is used to read data and perform
# initial filtering. It is also used to
# produce DTMs from input data using its GDAL
# driver
import pdal

# GDAL is used to write geoTIFFs for
# each product. OSR does some CRS handling
from osgeo import gdal
import osgeo.osr as osr

# geopandas gives us convenient data access
# and a spatial index
import geopandas as gpd
#import pandas as pd

# osmnx is used to conveniently create a fishnet
# grid of geometries at output resolution. It could
# be replaced by pure shapely if dependencies are an
# issue
import osmnx as ox

import os

# not using this, using geopandas instead
#from rtree import index


# two global variables....
THRESHOLD_HEIGHTS =[0.05, 0.10, 0.25, 0.5, 1, 2, 3]
NODATA_VALUE = -9999

# Vegetation cover fraction: (Nfirst - Nsingle) / Nfirst
def comp_vcf(points):
    """
    Computes vegetation cover fraction according to the TERN product manual.

    inputs:
    - a labelled array of points from an input LAS tile

    outputs:
    - a numpy array of grid cells containing the result of:

    (Nfirst - Nsingle) / Nfirst

    ...where:
    Nfirst = count of first returns
    Nsingle = count of single returns

    ...per grid cell.
    """
    # collect all the first and single return indices
    nSingle = np.size(np.where(points["NumberOfReturns"].values == 1))
    nFirst = np.size(np.where(points["ReturnNumber"].values == 1))
    if (nFirst > 0):
        vcf = (nFirst - nSingle) / nFirst
    else:
        print('no first returns, set vcf to {}'.format(NODATA_VALUE))
        vcf = NODATA_VALUE

    return(vcf)

    # Canopy layering index:

# (sum(weight * MaxNumberOfreturns * numberofpointsatreturn) / sum(weight * MaxNumberofReturns)) - 1

def comp_cli(points):
    maxreturns = np.max(points["NumberOfReturns"])
    minreturns = np.min(points["NumberOfReturns"])

    cli_numerator = []
    cli_denominator = []

    for nreturn in np.arange(minreturns, maxreturns+1, 1):
        # I don't really understand the fixed values in the weighting function - we can record
        # a lot more returns... I get that it's supposed to give higher return numbers a lower
        # weight based on likelihood of being recorded..

        weight = 1/min([5,nreturn])

        points_at_return = np.where(points["ReturnNumber"] == nreturn)
        cli_numerator.append(weight * maxreturns * np.size(points_at_return) )
        cli_denominator.append( weight * np.size(points_at_return ))

    cli = (np.sum(cli_numerator) / np.sum(cli_denominator)) - 1

    return(cli)

# vegetation layer cover fraction: LCF

# step 1: compute veg_below layers...

def comp_veg_layers(points, heights):
    """
    Compute veg layers needed for LCF as per the TERN product manual:

    LCF = VCF * (((veg returns below H2) - (veg returns below H1)) / (veg returns below H2))

    Inputs:
    - a set of points to compute vegetation counts over
    - an array of heights

    Outputs:
    - a dictionary of counts for vegetation points below each height

    Conditions:

    The LCF *must* be computed over the same set of points as the VCF used as input.

    """
    veg_below = {}
    #find veg returns - ASPRS classes 3,4,5
    veg_returns = np.where(np.logical_or(points["Classification"].values == 3,
                             points["Classification"].values == 4,
                             points["Classification"].values == 5))
    # set total of veg returns
    veg_below["all"] = np.size(veg_returns)

    # for each height threshold
    for height in heights:
        try:
            #get vegetation points in this cell below a given height threshold
            veg_below_height = np.size(np.where(points["HeightAboveGround"].values[veg_returns] < height))
        except ValueError:
            veg_below_height = np.nan
        # add the resulting value to a dictionary of point counts in this cell
        veg_below[str(height)] = veg_below_height

    return(veg_below)

def comp_lcf(veg_below, vcf):
    """
    Compute LCF as per the TERN product manual:

    LCF = VCF * (((veg returns below H2) - (veg returns below H1)) / (veg returns below H2))

    Inputs:
    - a dictionary of arrays with point counts per cell of vegetation below
     threshold heights, computed by comp_veg_below
    - a precomputed VCF raster

    Outputs:
    - a dictionary of arrays containing TERN manual layer fraction metrics

    Note: this function requires output of comp_veg_layers to be compiled into grids, it does not
          operate directly on the output of comp_veg_layers
    """
    lcf = {}

    #fuel layers
    lcf["lcf_cf"] = vcf * np.divide((veg_below["all"] - veg_below["2"]),
                        veg_below["all"])

    lcf["lcf_ef"] = vcf * np.divide((veg_below["2"] - veg_below["0.5"]),
                        veg_below["2"])

    lcf["lcf_nsf"] = vcf * np.divide((veg_below["0.5"] - veg_below["0.05"]),
                        veg_below["0.5"])

    #under and overstory
    lcf["lcf_h"] = vcf * np.divide((veg_below["1"] - veg_below["0.05"]),
                        veg_below["1"])
    lcf["lcf_os"] = vcf * np.divide((veg_below["3"] - veg_below["1"]),
                        veg_below["3"])
    lcf["lcf_us"] = vcf * np.divide((veg_below["all"] - veg_below["3"]),
                        veg_below["all"])

    return(lcf)

#CTH - TERN product manual says 'highest vegetation point', MATLAB code says '0.95 quantile of
# vegetation returns above 2m'. Below 2m is ignored

def comp_cth(points):
    # compute the highest vegetation point in each grid cell

    veg_returns = np.where(np.logical_or(points["Classification"].values == 3,
                             points["Classification"].values == 4,
                             points["Classification"].values == 5))

    try:

        vegpoints = points["HeightAboveGround"].values[veg_returns]
        canopy_index = np.where(vegpoints > 2.0)

        if (np.size(canopy_index) > 0):
            cth = np.quantile(vegpoints[canopy_index], 0.95)
        else:
            cth = NODATA_VALUE

    except ValueError:
        #print('no vegetation returns were present, CTH set to {}'.format(NODATA_VALUE))
        cth = NODATA_VALUE

    return(cth)

#CBH - ambiguous in TERN docs, will pull from MATLAB code
# there, it states that CBH is the 0.1th percentile of vegetation with normalised height
# above 2m
def comp_cbh(points):
    # compute the canopy base height in each cell.

    # grab an index of vegetation returns
    veg_returns = np.where(np.logical_or(points["Classification"].values == 3,
                             points["Classification"].values == 4,
                             points["Classification"].values == 5))
    try:
        #create an array of vegetation point normalised heights
        vegpoints = points["HeightAboveGround"].values[veg_returns]
        canopy_index = np.where(vegpoints > 2.0)
        if(np.size(canopy_index) > 0 ):
            #find the 0.1 quantile
            cbh = np.quantile(vegpoints[canopy_index], 0.10)
        else:
            cbh = NODATA_VALUE

    except ValueError:
        #if there are no veg points, set cth to NODATA

        #print('no vegetation returns were present, CTH set to {}'.format(NODATA_VALUE))
        cbh = NODATA_VALUE

    return(cbh)

def comp_fbf(points):
    # if building classes exist, compute a fractional cover per grid cell...
    # if no buildings exist return 0

    building_returns = np.where(points["Classification"].values == 6)

    totalpoints = np.size(points["Classification"])

    buildingpoints = np.size(building_returns)

    if (buildingpoints > 0):
        fbf = buildingpoints/totalpoints
    else:
        fbf = 0

    return(fbf)

def comp_density(points, resolution):
    """
    compute mean point density - npoints / area
    """

    #find number of points, any dimension will work here...
    npoints = np.size(points["Classification"])

    # divide npoints by area
    mean_density = npoints / resolution**2

    return(mean_density)

def comp_vh(points):
    """
    find vegetation points, compute mean height above ground for any vegetation
    present

    """

    try:
        veg_returns = np.where(np.logical_or(points["Classification"].values == 3,
                                 points["Classification"].values == 4,
                                 points["Classification"].values == 5))
        vh = np.mean(points["HeightAboveGround"].values[veg_returns])
    except ValueError:
        vh = np.nan

    return(vh)
