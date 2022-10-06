import os
from osgeo import gdal
from temp_path import splitPath, deleteFile, tempPath
from raster_options import RasterOptions, checkBandsMatch
from vector import VectorClipper


def warpRaster(rpath, tpath, **par):

    par = RasterOptions(**par)
    rds = gdal.Open(rpath)

    if rds:

        tsrs = par.getSpatialReference(rpath=rpath)
        if tsrs is None:
            tsrs = rds.GetSpatialRef()
        par['SpatialReference'] = tsrs

        if par.get('__preserve_original_pixel_size') and \
                (par.gget('PixelSize', 'PixelXSize', 'PixelYSize') is None):
            srs = rds.GetSpatialRef()
            if srs.IsProjected() and tsrs.IsProjected():
                trans = rds.GetGeoTransform()
                srs_unit_coef = srs.GetLinearUnits() / tsrs.GetLinearUnits()
                par.update(PixelXSize = abs(trans[1]*srs_unit_coef),
                           PixelYSize = abs(trans[5]*srs_unit_coef))

        if par.get('__vector_clipper'):
            vector_clipper = par['__vector_clipper']
            assert isinstance(vector_clipper, VectorClipper)
            par.update(vector_clipper.clipRasterParameters(rpath))

    rds = None

    return par.GDALWarpRaster(rpath, tpath)


def translateRaster(rpath, tpath, **parameters):
    return RasterOptions(**parameters).GDALTranslateRaster(rpath, tpath)


def rasterizeVector(vpath, tpath, **par):
    return RasterOptions(**par).GDALRasterize(vpath, tpath)


def saveRaster(rpath, tpath, **options):

    assert rpath != tpath
    overwrite = options.pop('overwrite')

    if os.path.exists(tpath):
        if overwrite:
            deleteFile(tpath)
        else:
            raise Exception('File already exists: ' + tpath)

    if options.get('__warp'):
        if checkBandsMatch(gdal.Open(rpath), options.get('Bands')) and \
                (not any([(key in options) for key in
                          ['mask', 'scale', 'exponent', 'colorinterp']])):
            warpRaster(rpath, tpath, **options)
        else:
            temp_path = tempPath(ext='tif', name=splitPath(rpath)[1])
            warp_options = RasterOptions(**options).alter(
                DataType=None, COMPRESS=None, NUM_THREADS='ALL_CPUS')
            warpRaster(rpath, temp_path, **warp_options)
            translateRaster(temp_path, tpath, **options)
            deleteFile(temp_path)
    else:
        translateRaster(rpath, tpath, **options)

