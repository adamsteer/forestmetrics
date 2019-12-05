#! /usr/bin/python3
"""
Utility functions for generating TERN forest products
from classified lidar point clouds

Dr Adam Steer
November 2019
"""

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
