# gdal_format
Python TUI for performing operations with raster files using GDAL functions with minimum code writing. Might make Python scripts and QGIS plugins a little shorter.

It is inspired with my 3 years of working as a GIS Developer in Python when I faced a lot of routine operations similar to each other. I needed to write the same functions (with the same typos) over and over again. So I made several functions to make my life easier. Here I tried to put 'em all together hoping it would be useful for anybody else.

Designed for making GeoTIFF, but can use most of file types GDAL can open.

Most of GDAL options can be written into raster_options.RasterOptions() class like this:

    options = RasterOptions(COMPRESS='DEFLATE', ZLEVEL=9, PREDICTOR=2, NUM_THREADS='ALL_CPUS')

Some raster parameters may be set through new options like:

    options = RasterOptions(PixelSize=100, DataType=1, NoDataValue=0)

Also some special options are used, like:

__preserve_original_pixel_size - if True pixel size is not changing when transformed between two projective CRS

__wrap - force using gdal_warp in saveRaster()

__vector_clipper - used for clipping raster files with vector to make sure the target raster extent corresponds to the intersection between the source raster and vector. Must contain vector.VectorClipper(<path_to_the_cutline>) object. If raster is within vector no crop operation is done.

The key functions to manipulate raster are in save_data.py file and they are the following:

    warpRaster(rpath, tpath, **options) # warp_raster with the new options
  
    translateRaster(rpath, tpath, **options) # translate_raster with the new options
  
    saveRaster(rpath, tpath, **options) # uses translate_raster and warp_raster to make the very raster file you want (make sure __warp=True to use gdal_warp)

Example:
   
    from raster_options import RasterOptions
    from save_data import saveRaster
    
    path_in = r'c:\test\raster_in.tif'
    path_out = r'c:\test\raster_out.tif'
    path_clip = r'c:\test\clip.shp'
    
    options = RasterOptions(EPSG=3857, PixelSize=100, Bands=[3,2,1], NUM_THREADS='ALL_CPUS', __warp=True, overwrite=True)
    options.makeVectorClipper(vector_path=path_clip)
    
    saveRaster(path_in, path_out, **options)

Creates a new raster file with 3 bands corresponding to bands 3,2,1 of the original raster, reprojected to EPSG:3857, with PixelSize=100m, and clipped with r'c:\test\clip.shp'
