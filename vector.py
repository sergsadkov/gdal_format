from osgeo import gdal, ogr, osr
from temp_path import TempPathsStorage, tempPath


# Returns geometry fully covering vector layer geometry
def vectorLayerGeometry(data_source, iLayer=0, geom_type=ogr.wkbMultiPolygon):
    layer = data_source.GetLayer(iLayer)
    if layer is not None:
        vector_geometry = ogr.Geometry(geom_type)
        for feature in layer:
            feature_geometry = feature.GetGeometryRef()
            if feature_geometry.GetGeometryType() == ogr.wkbMultiPolygon:
                vector_geometry = vector_geometry.Union(feature_geometry)
            else:
                vector_geometry.AddGeometry(feature_geometry)
        return vector_geometry


def rasterDataCover(rpath, vpath, mask_band_id=1):

    # Mapping between gdal type and ogr field type
    type_mapping = {gdal.GDT_Byte: ogr.OFTInteger,
                    gdal.GDT_UInt16: ogr.OFTInteger,
                    gdal.GDT_Int16: ogr.OFTInteger,
                    gdal.GDT_UInt32: ogr.OFTInteger,
                    gdal.GDT_Int32: ogr.OFTInteger,
                    gdal.GDT_Float32: ogr.OFTReal,
                    gdal.GDT_Float64: ogr.OFTReal,
                    gdal.GDT_CInt16: ogr.OFTInteger,
                    gdal.GDT_CInt32: ogr.OFTInteger,
                    gdal.GDT_CFloat32: ogr.OFTReal,
                    gdal.GDT_CFloat64: ogr.OFTReal}

    ds = gdal.Open(rpath)
    prj = ds.GetProjection()
    src_band = ds.GetRasterBand(mask_band_id)
    mask_band = src_band.GetMaskBand()

    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(vpath)
    srs = osr.SpatialReference(wkt=prj)
    dst_layer = dst_ds.CreateLayer('', srs=srs)

    raster_field = ogr.FieldDefn('id', type_mapping[src_band.DataType])
    dst_layer.CreateField(raster_field)
    options = ['DATASET_FOR_GEOREF=' + rpath]

    gdal.Polygonize(mask_band, mask_band, dst_layer, 0, options, callback=None)

    dst_ds = None


def rasterCoverGeometry(raster_path, raster_cover_path=None, srs=None):

    temp_paths = TempPathsStorage()

    if raster_cover_path is None:
        raster_cover_path = temp_paths.newTempPath(ext='shp')

    rasterDataCover(raster_path, raster_cover_path)
    cover_geometry = vectorLayerGeometry(ogr.Open(raster_cover_path),
                                         geom_type=ogr.wkbMultiPolygon)

    if srs is not None:
        cover_geometry = geometryCoordinateTransformation(
            cover_geometry, gdal.Open(raster_path).GetSpatialRef(), srs)

    del temp_paths

    return cover_geometry


def geometryCoordinateTransformation(geom, start_srs, end_srs):
    if (start_srs is not None) and (end_srs is not None):
        if start_srs.ExportToProj4() != end_srs.ExportToProj4():
            geom.Transform(osr.CoordinateTransformation(start_srs, end_srs))
    return geom


def geometryExtent(geometry):
    xmin, xmax, ymin, ymax = geometry.GetEnvelope()
    return [xmin, ymin, xmax, ymax]


def geometryShapefile(geometry, vector_path, srs=None):

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dsc = driver.CreateDataSource(vector_path)

    if dsc is not None:
        layer = dsc.CreateLayer('', srs=srs, geom_type=geometry.GetGeometryType())
        feature = ogr.Feature(ogr.FeatureDefn())
        feature.SetGeometry(geometry)
        layer.CreateFeature(feature)
        dsc = None

        prj_path = vector_path[:-4] + '.prj'
        with open(prj_path, 'w') as prj:
            prj.write(srs.ExportToWkt())


# Contains data for clipping raster with vector
class VectorClipper(object):

    def __init__(self, vector_path):

        self.path = vector_path
        self.dsc = ogr.Open(vector_path)

        if self.dsc is None:
            self.layer = None
            self.srs = None
            self.geometry = None
        else:
            self.layer = self.dsc.GetLayer()
            self.srs = self.layer.GetSpatialRef()
            self.geometry = vectorLayerGeometry(self.dsc,
                                                geom_type=ogr.wkbMultiPolygon)

    def __str__(self):
        if self.srs is None:
            return f'VectorClipper: "{self.path}" without spatial reference data'
        else:
            return f'VectorClipper: {self.path} in {self.srs.ExportToProj4()}'


    def __repr__(self):
        return self.__str__()


    def intersect(self, check_geometry):
        if self.geometry.Intersects(check_geometry):
            if check_geometry.Within(self.geometry):
                return 2
            elif self.geometry.Within(check_geometry):
                return 3
            else:
                return 1
        else:
            return 0


    def intersection(self, check_geometry):
        intersect_status = self.intersect(check_geometry)
        if intersect_status == 0:
            return None
        elif intersect_status == 1:
            return self.geometry.Intersection(check_geometry)
        elif intersect_status == 2:
            return check_geometry
        elif intersect_status == 3:
            return self.geometry
        else:
            raise Exception(f'Wrong intersect status {intersect_status}, a value between 0 and 3 is required')


    def clipRasterParameters(self, raster_path, end_srs=None):

        temp_paths = TempPathsStorage()

        raster_cover_path = temp_paths.newTempPath(ext='shp', name='tmp')
        rasterDataCover(raster_path, raster_cover_path)
        cover_geometry = vectorLayerGeometry(ogr.Open(raster_cover_path),
                                             geom_type=ogr.wkbMultiPolygon)
        cover_geometry = geometryCoordinateTransformation(
            cover_geometry, gdal.Open(raster_path).GetSpatialRef(), self.srs)

        intersect_status = self.intersect(cover_geometry)

        if intersect_status == 0:
            # No intesection -- no raster saved
            result = f'{self.path} and {raster_path} areas do not intersect, cannot clip'

        elif intersect_status == 1:
            # Raster intersects vector -- intersection vector is used
            intersection_geometry = self.geometry.Intersection(cover_geometry)
            if end_srs is not None:
                intersection_geometry = geometryCoordinateTransformation(
                    intersection_geometry, self.srs, end_srs)
            clip_geometry_path = tempPath(ext='shp', name='tmp')
            print(clip_geometry_path, intersection_geometry.ExportToWkt()[:100])
            geometryShapefile(intersection_geometry, clip_geometry_path,
                              srs=self.srs)
            import os
            print(os.path.exists(clip_geometry_path))

            result = {'crop_to_cutline': True,
                      'cutline': clip_geometry_path,
                      'Extent': geometryExtent(intersection_geometry)}

        elif intersect_status == 2:
            # Raster within vector -- no raster cut
            result = {'crop_to_cutline': False,
                      'cutline': None,
                      'Extent': None}

        elif intersect_status == 3:
            # Vector within raster -- the original cut vector is used
            result = {'crop_to_cutline': True,
                      'cutline': self.path,
                      'Extent': geometryExtent(self.geometry)}

        # Wrong intersect status
        else:
            result = f'Wrong intersect status {intersect_status}, a value between 0 and 3 is required'

        # del temp_paths

        if isinstance(result, str):
            raise Exception(result)
        else:
            return result
