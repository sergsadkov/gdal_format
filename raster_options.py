import re
import math
from copy import copy
from osgeo import gdal, ogr, osr


# Does not support the following gdal_translate options:
#   -tr -expand -outsize -srcwin -projwin -projwin_srs -a_srs -a_coord_epoch
#   -a_offset -a_ullr -colorinterp_X -gcp -sds -norat -noxmp

# Does not support the following gdal_rasterize options:
#   -of -a_srs -te -tap -ts


__default = {
        'RasterSize_error': 0.01,
        'PixelSize_error': 0.01,
        'NoDataValue_error': 0.01,
        'DataMinimum_error': 0.01,
        'DataMaximum_error': 0.01,
        '__preserve_original_pixel_size': True,
        '__vector_clipper': None,
        '__wrap': False,
    }

__gdal_data_type = ['Unknown', 'Byte', 'UInt16', 'Int16', 'UInt32', 'Int32',
                    'Float32', 'Float64', 'CInt32', 'CFloat32', 'CFloat64']

__int_options = ['RasterCount', 'RasterCountMin', 'RasterCountMax',
                    'RasterSize', 'RasterSizeMin', 'RasterSizeMax',
                    'DataType', 'EPSG']

__float_options = ['PixelSize', 'PixelSizeMin', 'PixelSizeMax',
                      'NoDataValue', 'NoDataValueMin', 'NoDataValueMax',
                      'DataMinimum', 'DataMinimumMin', 'DataMinimumMax',
                      'DataMaximum', 'DataMaximumMin', 'DataMaximumMax']

__string_options = ['Projection', 'Compression']

# See details on https://gdal.org/drivers/raster/gtiff.html#creation-options
__gdal_co = ['NUM_THREADS', 'GEOREF_SOURCES', 'SPARSE_OK', 'TWF', 'RPB',
             'RPCTXT', 'INTERLEAVE', 'TILED', 'BLOCKXSIZE', 'BLOCKYSIZE',
             'NBITS', 'COMPRESS', 'PREDICTOR', 'DISCARD_LSB', 'JPEG_QUALITY',
             'JPEGTABLESMODE', 'ZLEVEL', 'ZSTD_LEVEL', 'MAX_Z_ERROR',
             'WEBP_LEVEL', 'WEBP_LOSSLESS', 'JXL_LOSSLESS', 'JXL_EFFORT',
             'JXL_DISTANCE', 'PHOTOMETRIC', 'ALPHA', 'PROFILE', 'BIGTIFF',
             'PIXELTYPE', 'COPY_SRC_OVERVIEWS', 'GEOTIFF_KEYS_FLAVOR',
             'GEOTIFF_VERSION']

# See details on https://gdal.org/api/gdalwarp_cpp.html#_CPPv4N15GDALWarpOptions16papszWarpOptionsE
__gdal_wo = ['INIT_DEST', 'WRITE_FLUSH', 'SKIP_NOSOURCE', 'UNIFIED_SRC_NODATA',
             # 'CUTLINE',
             'CUTLINE_BLEND_DIST', 'CUTLINE_ALL_TOUCHED', 'OPTIMIZE_SIZE',
             'NUM_THREADS', 'STREAMABLE_OUTPUT', 'SRC_COORD_PRECISION',
             'SRC_ALPHA_MAX', 'DST_ALPHA_MAX', 'SAMPLE_GRID', 'SAMPLE_STEPS',
             'SOURCE_EXTRA', 'APPLY_VERTICAL_SHIFT',
             'MULT_FACTOR_VERTICAL_SHIFT']

# See details on https://gdal.org/api/gdal_alg.html#_CPPv426GDALCreateRPCTransformerV2PK13GDALRPCInfoV2idPPc
__gdal_to = ['RPC_HEIGHT', 'RPC_HEIGHT_SCALE', 'RPC_DEM',
             'RPC_DEM_INTERPOLATION', 'RPC_DEM_MISSING_VALUE', 'RPC_DEM_SRS',
             'RPC_DEM_APPLY_VDATUM_SHIFT', 'RPC_PIXEL_ERROR_THRESHOLD',
             'RPC_MAX_ITERATIONS', 'RPC_FOOTPRINT']

# See details on https://gdal.org/drivers/raster/gtiff.html
__gdal_oo = ['NUM_THREADS', 'GEOREF_SOURCES', 'SPARSE_OK']


__gdal_options = {'creation': __gdal_co, 'transform': __gdal_to,
                'warp': __gdal_wo, 'open': __gdal_oo}


class GDALoptions:

    def __init__(self, gdalf):

        if gdalf == 'gdalwarp':
            self.single = ['tps', 'rpc', 'geoloc', 'tap', 'multi', 'q',
                           'overwrite', 'crop_to_cutline', 'nomd', 'setci',
                           'srcalpha', 'nosrcalpha', 'dstalpha']
            self.string = ['ovr', 'if', 'of', 'cutline', 'cl', 'cwhere',
                           'csql', 'cvmd']
            self.integer = ['order', 'refine_gcps', 'wm', 'cblend']
            self.float = ['et', 'srcnodata', 'dstnodata']

        elif gdalf == 'gdal_translate':
            self.single = ['strict', 'unscale', 'epo', 'eco', 'nogcp', 'sds',
                           'q', 'stats', 'norat', 'noxmp', 'overwrite']
            self.string = ['if', 'of', 'expand', 'colorinterp', 'scale']
            self.integer = ['mask']
            self.float = ['exponent', 'a_scale', 'a_offset', 'a_nodata']

        elif gdalf == 'gdal_rasterize':
            self.single_options = ['i', 'at', '3d', 'add', 'q', 'overwrite']
            self.str_options = ['a', 'l', 'sql', 'dialect', 'init', 'optim']
            self.int_options = ['burn']
            self.float_options = []

        else:
            raise Exception(f'Unknown GDAL function: {gdalf}')


def toList(obj):
    new_obj = copy(obj)
    if isinstance(new_obj, (tuple, list)):
        return new_obj
    else:
        return [new_obj]


def checkVar(var, val_exact, val_presicion, val_min, val_max):

    if val_exact is not None:

        if hasattr(var, '__float__') and (val_presicion is not None):

            try:
                if abs((var - val_exact) / (var + val_exact)) > val_presicion:
                    return 1

            except ZeroDivisionError:
                if abs(var) > val_presicion:
                    return 1

        elif var != val_exact:
            return 1

    if hasattr(var, '__float__'):

        if val_min is not None:
            if var < val_min:
                return 2

        if val_max is not None:
            if var > val_max:
                return 3

    return 0


def getSpatialReferenceFromEPSG(epsg):
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))
    return srs


def getSpatialReferenceFromProj4(proj4):
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj4)
    return srs


def getUTMzoneEPSG(lat, lon):
    return int(f'32{(6,7)[lat<0]}{round(((lon+180)%360)/6):02}')


def getRasterGeometry(rds):

    rgt = rds.GetGeoTransform(can_return_null=True)

    if rgt is not None:
        lx = rgt[0]
        ly = rgt[3]
        hx = rgt[0] + (rgt[1] * math.cos(rgt[2]) * rds.RasterXSize) + \
             (rgt[5] * math.sin(rgt[4]) * rds.RasterYSize)
        hy = rgt[3] + (rgt[5] * math.cos(rgt[4]) * rds.RasterYSize) + \
             (rgt[1] * math.sin(rgt[2]) * rds.RasterXSize)
        wkt = f'POLYGON (({hx} {hy}, {hx} {ly}, {lx} {ly}, {lx} {hy}, {hx} {hy}))'

    else:
        point_coords = ["%f %f" % (p.GCPX, p.GCPY) for p in rds.GetGCPs()]
        wkt = f'MULTIPOINT ({",".join(point_coords)})'

    return ogr.Geometry(wkt=wkt)


def getUTMSpatialReference(rds):

    rsrs = rds.GetSpatialRef()

    if rsrs is None:
        raise Exception('Spatial reference not found, cannot find UTM zone')

    tsrs = getSpatialReferenceFromEPSG(4326)

    if rsrs.ExportToProj4() == tsrs.ExportToProj4():
        transform = None
    else:
        # Need to AxisMappingStrategy for different sources
        rsrs.SetAxisMappingStrategy(0)
        tsrs.SetAxisMappingStrategy(0)
        transform = osr.CoordinateTransformation(rsrs, tsrs)

    area_geometry = getRasterGeometry(rds)
    if transform:
        area_geometry.Transform(transform)
    point_geometry = area_geometry.Centroid()
    utm_epsg = getUTMzoneEPSG(point_geometry.GetY(), point_geometry.GetX())

    return getSpatialReferenceFromEPSG(utm_epsg)


def checkBandsMatch(rds, bands):
    if bands is None:
        return True
    if rds is not None:
        return list(range(1, rds.RasterCount+1)) == bands
    else:
        raise Exception("Raster dataset not found, cannot check")


def unbanded(option):
    if re.search('^.+_\d+$', option):
        return option.split('_')[0]
    else:
        return option


def getPercentile(histogram, min=0.02, max=0.98):
    total = sum(histogram)
    min_value = total * min
    max_value = total * max
    min_position = 0
    min_sum = histogram[0]
    while min_sum < min_value:
        min_position += 1
        min_sum += histogram[min_position]
    max_position = len(histogram) - 1
    max_sum = total - histogram[max_position]
    while max_sum > max_value:
        max_position -= 1
        max_sum -= histogram[max_position]
    return min_position, max_position


class RasterOptions(dict):

    def __init__(self, *args, **kwargs):
        self.update(globals()['__default'])
        for arg in args:
            self.update(arg)
        self.update(**kwargs)


    def gget(self, *keys, default=None, return_key=False):

        for key in keys:
            if key in self:
                if return_key:
                    return key, self[key]
                else:
                    return self[key]

        if return_key:
            return None, default
        else:
            return default


    def valueList(self, *keys, default=None):
        return [self.get(key, default) for key in keys]


    def alter(self, **kwargs):
        new_dict = copy(self)
        new_dict.update(**kwargs)
        return new_dict


    def getSpatialReference(self, rpath=None, for_extent=False):

        dkey = ('', 'Extent')[for_extent]
        srs_keys = ['SpatialReference', 'Proj4', 'ProjWkt', 'EPSG']
        keys = tuple([dkey + key for key in srs_keys])
        key, srs = self.gget(*keys, return_key=True)

        if (key is not None) and (srs is not None):

            if srs == 'UTM':
                if rpath is not None:
                    return getUTMSpatialReference(gdal.Open(rpath))
                else:
                    return srs

            elif key == 'SpatialReference':
                return srs

            elif key == 'Proj4':
                return getSpatialReferenceFromProj4(srs)

            elif key == 'ProjWkt':
                return osr.SpatialReference(wkt = srs)

            elif key == 'EPSG':
                return getSpatialReferenceFromEPSG(srs)

            else:
                raise Exception('Unknown srs format:', srs)


    def getGDALoptions(self, gdalf):

        options = []
        gdalo = GDALoptions(gdalf)

        for option in self:

            value = self.get(option)

            if value is not None:

                if unbanded(option) in gdalo.single:
                    if value:
                        options.append(f' -{option}')

                elif unbanded(option) in gdalo.string:
                    options.append(f' -{option} {value}')

                elif unbanded(option) in (gdalo.integer + gdalo.float):
                    options.append(f' -{option} {value}')

        return options


    def getGDALoptionsList(self, type):

        options = []

        gdal_options = globals()['__gdal_options'].get(type, [])

        for option in self:
            option_check = option.upper()
            if option_check in gdal_options:
                options.append(f'{option_check}={str(self[option])}')

        return options


    # Possible gdalf values: gdalwarp, gdal_translate
    def getGDALoptionString(self, gdalf):

        options = self.getGDALoptions(gdalf)

        if self.get('DataType') is not None:
            options.append(f' -ot {__gdal_data_type[self.get("DataType", 0)]}')

        if self.get('Method') is not None:
            options.append(f' -r {self["Method"]}')

        if self.get('NoDataValue') is not None:
            nodata_apx = {"gdalwarp": "dst",
                          "gdal_translate": "a_",
                          "gdal_rasterize": "a_"}.get(gdalf, "")
            options.append(f' -{nodata_apx}nodata {self["NoDataValue"]}')

        if gdalf == 'gdalwarp':

            srs = self.getSpatialReference()

            if srs is not None:
                if isinstance(srs, osr.SpatialReference):
                    options.append(f' -t_srs "{srs.ExportToProj4()}"')

            if self.get('Extent') is not None:
                extent_str = " ".join([str(i) for i in self["Extent"]])
                options.append(f' -te {extent_str}')
                extent_srs = self.getSpatialReference(for_extent=True)
                if isinstance(extent_srs, osr.SpatialReference):
                    options.append(f' -te_srs "{extent_srs.ExportToProj4()}"')

            rasterXsize = self.gget('RasterXSize', 'RasterSize')
            rasterYsize = self.gget('RasterYSize', 'RasterSize')

            if rasterXsize and rasterYsize:
                options.append(f' -ts {rasterXsize} {rasterYsize}')

            pixelXsize = self.gget('PixelXSize', 'PixelSize')
            pixelYsize = self.gget('PixelYSize', 'PixelSize')

            if pixelXsize and pixelYsize:
                options.append(f' -tr {pixelXsize} {pixelYsize}')

            if self.get('SourceNoDataValue') is not None:
                options.append(f' -srcnodata {self["SourceNoDataValue"]}')

        elif gdalf == 'gdal_translate':

            if self.get('Bands') is not None:
                options.append(''.join([f' -b {int(bandnum)}'
                                        for bandnum in self.get('Bands')]))

            if self.get('mask') is not None:
                options.append(' --config GDAL_TIFF_INTERNAL_MASK YES')

        elif gdalf == 'gdal_rasterize':
            if sum([key in self for key in ['3d', 'burn', 'a']]) != 1:
                raise Exception('One and only one of "3d", "burn" or "a" is required')

        option_type = {'gdalwarp': ['creation', 'transform', 'warp', 'open'],
                       'gdal_translate': ['creation', 'open'],
                       'gdal_rasterize': ['creation', 'open'],
                       }.get(gdalf, [])

        for type in option_type:
            options += [f' -{type[0].lower()}o {option}'
                        for option in self.getGDALoptionsList(type)]

        return ''.join(options)


    def GDALTranslateRaster(self, rpath, tpath):

        rds = gdal.Open(rpath)
        option_string = self.getGDALoptionString('gdal_translate')

        # if option GDAL_TIFF_INTERNAL_MASK set to YES
        if self.get('mask') is not None:
            option_string = option_string.replace(
                ' --config GDAL_TIFF_INTERNAL_MASK YES', '')
            gdal.SetConfigOption('GDAL_TIFF_INTERNAL_MASK', 'YES')

        translate_options = gdal.TranslateOptions(gdal.ParseCommandLine(option_string))
        gdal.Translate(tpath, rds, options=translate_options)
        rds = None


    def GDALWarpRaster(self, rpath, tpath):
        option_string = self.getGDALoptionString('gdalwarp')
        warp_options = gdal.WarpOptions(gdal.ParseCommandLine(option_string))
        gdal.Warp(tpath, rpath, options=warp_options)


    def GDALRasterize(self, vpath, tpath):
        option_string = self.getGDALoptionString('gdal_rasterize')
        rasterize_options = gdal.RasterizeOptions(gdal.ParseCommandLine(option_string))
        gdal.Rasterize(gdal.Open(tpath, 1), vpath, options=rasterize_options)
